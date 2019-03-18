#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

double ceiled(double f)
{
    return ceil(f-0.00001);
}

double floored(double f)
{
    return floor(f+0.00001);
}


vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      unsigned char color[3];
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;
};

std::vector<Triangle>
GetTriangles(void)
{
   std::vector<Triangle> rv(100);

   unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
   for (int i = 0 ; i < 100 ; i++)
   {
       int idxI = i%10; // 0,1,2,3,4,5,6,7,8,9
       int posI = idxI*100; //0,100,200,300,400,500,600,700,800,900
       int idxJ = i/10;  //ten 0s, ten 1s, ten 2s, ten 3s, ten 4s,ten 5s,ten 6s,ten 7s,ten 8s,ten 9s
       int posJ = idxJ*100; //ten 0s, ten 100s, ten 200s, ten 300s, ten 400s,ten 500s,ten 600s,ten 700s,ten 800s,ten 900s
       int firstPt = (i%3); //0,1,2,0,1,2,...,0,1,2,0
       rv[i].X[firstPt] = posI;  //the X value of ith tirangle's 1st vertex
       if (i == 50)
           rv[i].X[firstPt] = -10;  //the X value of 50th tirangle's 1st vertex is -10
       rv[i].Y[firstPt] = posJ+10*(idxJ+1); //the Y value of ith tirangle's 1st vertex
       rv[i].X[(firstPt+1)%3] = posI+99;  //the X value of ith tirangle's 2nd vertex
       rv[i].Y[(firstPt+1)%3] = posJ+10*(idxJ+1);  //the Y value of ith tirangle's 2nd vertex
       rv[i].X[(firstPt+2)%3] = posI+i;  //the X value of ith tirangle's 3rd vertex
       rv[i].Y[(firstPt+2)%3] = posJ;  //the Y value of ith tirangle's 3rd vertex
       if (i == 5)
          rv[i].Y[(firstPt+2)%3] = -50;	//the Y value of 5th tirangle's 3rd vertex is -50
       rv[i].color[0] = colors[i%6][0]; //the r channel for ith triangle
       rv[i].color[1] = colors[i%6][1]; //the g color channel for ith triangle
       rv[i].color[2] = colors[i%6][2]; //the b color channel for ith triangle
   }

   return rv;
}

double GetSlope(double x1, double x2, double y1, double y2) {
    if ((x1 - x2) != 0) {
        return (y1-y2)/(x1-x2);
    }
    else {
        return 0.00;
    }
}

double GetIntersectionValueX(double scanLine, double slope, double x, double y)
{
    if (slope != 0.00) {
        return (scanLine-y)/slope + x;
    }
    else 
        return x;
}

void RasterizeGoingDownTriangle(Triangle &t, Screen screen) {
    double x_left = t.X[0];
    double x_right = t.X[0];
    double x_down = t.X[0];
    double y_up = t.Y[0];
    double y_down = t.Y[0];

    for (int i=1; i<3; i++) {
        if (t.Y[i] > y_up) {
          y_up = t.Y[i];
        }
        if (t.Y[i] < y_down) {
          y_down = t.Y[i];
          x_down = t.X[i];
        }

        if (t.X[i] > x_right) {
          x_right = t.X[i];
        }
        if (t.X[i] < x_left) {
          x_left = t.X[i];
        }
    }

    double rowMin = ceiled(y_down);
    double rowMax = floored(y_up);
    if (rowMin <= rowMax) {
        for (int r = rowMin; r <= rowMax; r++) {
            double slopeLeft = GetSlope(x_left, x_down, y_up, y_down);
            double slopeRight = GetSlope(x_right, x_down, y_up, y_down);

            double leftEnd = GetIntersectionValueX(r, slopeLeft, x_down, y_down);
            double rightEnd = GetIntersectionValueX(r, slopeRight, x_down, y_down);

            leftEnd = ceiled(leftEnd);
            rightEnd = floored(rightEnd);

            for (int c = leftEnd; c <= rightEnd; c++) {
                if (r >= screen.height || c < 0 || c >= screen.width || r < 0) { continue; }
                else {
                    int index = r*screen.width + c;
                    screen.buffer[index*3] = t.color[0];
                    screen.buffer[index*3+1] = t.color[1];
                    screen.buffer[index*3+2] = t.color[2];                
                }
            }
        }        
    }
}
bool isGoingDown(Triangle &t) {
    double maxY = t.Y[0];
    for (int i = 0; i < 3; i++) {
        if (t.Y[i] > maxY) {
            maxY = t.Y[i];
        }
    }
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            if ((fabs(t.Y[i] - t.Y[j]) < 0.000000000001) && (fabs(t.Y[i] - maxY) < 0.000000000001)){
                return true;
            }
        }
    }
    return false;
}

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;

   // YOUR CODE GOES HERE TO DEPOSIT THE COLORS FROM TRIANGLES 
   // INTO PIXELS USING THE SCANLINE ALGORITHM
   for (int i = 0; i < triangles.size(); i++)
   {
    if (isGoingDown(triangles[i])) {
      RasterizeGoingDownTriangle(triangles[i], screen);
    }
   }
   WriteImage(image, "allTriangles");
}
