#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <tgmath.h>

using std::cerr;
using std::endl;

int triangleID = -1;

class Triangle
{
public:
    double              X[3];
    double              Y[3];
    unsigned char       color[3];
};

class Screen
{
public:
    double          width;
    double          height;
    unsigned char   *buffer;
    void fillColor(unsigned char *color, int c, int r) {
        if (r >= height || c < 0 || c >= width || r < 0) { return; }
        // cerr << "Triangle " << triangleID << " is writing to row " << r << ", column " << c << endl;
        int index = r*width+c;
        buffer[index*3] = color[0];
        buffer[index*3+1] = color[1];
        buffer[index*3+2] = color[2];
    }
};

int ceiled(double f)
{
    return ceil(f-0.00001);
}

int floored(double f)
{
    return floor(f+0.00001);
}

double GetSlope(double x1, double x2, double y1, double y2) {
    if (x1 != x2) {
        return (y1-y2)/(x1-x2);
    }
    else {
        return 99999999.00;
    }
}

double GetIntersectionValueX(int scanLine, double slope, double x, double y)
{
    double doublescanline = (double)scanLine;
    if (slope == 99999999.00) {
        return x;
    }
    else {
        return (doublescanline-y)/slope + x;
    }
}

vtkImageData *NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1c_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
    float *color_ptr = colors->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
        tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
        tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
    }
    cerr << "Done reading" << endl;

    return tris;
}

bool isArbitrary(Triangle &t) {
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            if (fabs(t.Y[i] - t.Y[j]) < 0.000000000001) {
                return false;
            }
        }
    }
    return true;
}

bool isGoingUp(Triangle &t) {
    double minY = t.Y[0];
    for (int i = 0; i < 3; i++) {
        if (t.Y[i] < minY) {
            minY = t.Y[i];
        } 
    }
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            if ((fabs(t.Y[i] - t.Y[j]) < 0.000000000001) && (fabs(t.Y[i] - minY) < 0.000000000001)){
                return true;
            }
        }
    }
    return false;
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

void RasterizeGoingDownTriangle(Triangle &t, Screen screen) {
    double x_left = 0.0;
    double x_right = 0.0;
    double x_down = t.X[0];
    double y_up_left = 0.0;
    double y_up_right = 0.0;
    double y_down = t.Y[0];
    
    for (int i=1; i<3; i++) {
        if (t.Y[i] < y_down) {
          y_down = t.Y[i];
          x_down = t.X[i];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            if (t.Y[i] == t.Y[j]) {
                if (t.X[j] >= t.X[i]) {
                    x_left = t.X[i];
                    y_up_left = t.Y[i];
                    x_right = t.X[j];
                    y_up_right = t.Y[j];
                }
                else {
                    x_left = t.X[j];
                    y_up_left = t.Y[j];
                    x_right = t.X[i];
                    y_up_right = t.Y[i];                
                }
            }
        }
    }

    int rowMin = ceiled(y_down);
    int rowMax = floored(y_up_left);
    for (int r = rowMin; r <= rowMax; r++) {
        double slopeLeft = GetSlope(x_left, x_down, y_up_left, y_down);
        double slopeRight = GetSlope(x_right, x_down, y_up_left, y_down);

        double leftEnd = GetIntersectionValueX(r, slopeLeft, x_down, y_down);
        double rightEnd = GetIntersectionValueX(r, slopeRight, x_down, y_down);

        int left = ceiled(leftEnd);
        int right = floored(rightEnd);
        for (int c = left; c <= right; c++) {
            screen.fillColor(t.color, c, r);
        }
    }
}

void RasterizeGoingUpTriangle(Triangle &t, Screen screen) {
    double x_up = t.X[0];
    double y_up = t.Y[0];
    double x_left = t.X[0];
    double x_right = t.X[0];
    double y_down_left = t.Y[0];
    double y_down_right = t.Y[0];

    for (int i=1; i<3; i++) {
        if (t.Y[i] > y_up) {
          y_up = t.Y[i];
          x_up = t.X[i];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            if (t.Y[i] == t.Y[j]) {
                if (t.X[j] >= t.X[i]) {
                    x_left = t.X[i];
                    y_down_left = t.Y[i];

                    x_right = t.X[j];
                    y_down_right = t.Y[j];
                }
                else {
                    x_left = t.X[j];
                    y_down_left = t.Y[j];

                    x_right = t.X[i];
                    y_down_right = t.Y[i];
                }
            }
        }
    }

    int rowMin = ceiled(y_down_right);
    int rowMax = floored(y_up);
    for (int r = rowMin; r <= rowMax; r++) {
        double slopeLeft = GetSlope(x_left, x_up, y_down_left, y_up);
        double slopeRight = GetSlope(x_right, x_up, y_down_right, y_up);

        double leftEnd = GetIntersectionValueX(r, slopeLeft, x_up, y_up);
        double rightEnd = GetIntersectionValueX(r, slopeRight, x_up, y_up);

        leftEnd = ceiled(leftEnd);
        rightEnd = floored(rightEnd);
        
        for (int c = leftEnd; c <= rightEnd; c++) {
            screen.fillColor(t.color, c, r);
        }
    }
}

void RasterizeArbitraryTriangle(Triangle &t, Screen screen)
{
    double y_down = t.Y[0];
    double y_up = t.Y[0];
    double x_down = t.X[0];
    double x_up = t.X[0];
    double x_mid = 0.0;
    double y_mid = 0.0;


    for (int i = 1; i < 3; i++) {
        if (t.Y[i] > y_up) { 
            y_up = t.Y[i];
            x_up = t.X[i];
        }
        if (t.Y[i] < y_down) {
            y_down = t.Y[i];
            x_down = t.X[i];
        }
    }

    for (int i = 0; i < 3; i++) {
        if ((t.Y[i] < y_up) && (t.Y[i] > y_down)) {
            y_mid = t.Y[i];
            x_mid = t.X[i];
        }
    }

    double slope_up_down = GetSlope(x_up, x_down, y_up, y_down);
    double slope_mid_down = GetSlope(x_mid, x_down, y_mid, y_down);
    double slope_mid_up = GetSlope(x_mid, x_up, y_mid, y_up);
    double x_mid_alt = (y_mid-y_up)/slope_up_down + x_up;

// Goling down
    int rowMinGoingDown = ceiled(y_down);
    int rowMaxGoingDown = floored(y_mid);
    for (int r = rowMinGoingDown; r <= rowMaxGoingDown; ++r) {
        double mid_down = GetIntersectionValueX( r, slope_mid_down, x_mid, y_mid);
        double up_down = GetIntersectionValueX(r, slope_up_down, x_mid_alt, y_mid);

        int leftEnd, rightEnd;
        if (mid_down <= up_down) {
            
            leftEnd = ceiled(mid_down);
            rightEnd = floored(up_down);
        }
        else {
            leftEnd = ceiled(up_down);
            rightEnd = floored(mid_down);
        }
        if (rightEnd - leftEnd > 20) {
            abort();
        }

        for (int c = leftEnd; c <= rightEnd; c++) {
            screen.fillColor(t.color, c, r);
        }
    }

// Going up
    int rowMinGoingUp = ceiled(y_mid);
    int rowMaxGoingUp = floored(y_up);
    int static count = 0;

    for (int row = rowMinGoingUp; row <= rowMaxGoingUp; ++row) {
        double up_mid = GetIntersectionValueX(row, slope_mid_up, x_up, y_up);
        double down_up = GetIntersectionValueX(row, slope_up_down, x_up, y_up);
        double left, right;
        if (up_mid >= down_up) {
            left = ceiled(down_up);
            right = floored(up_mid);
        }
        else {
            left = ceiled(up_mid);
            right = floored(down_up);
        }
        if (right - left > 20) {
            abort();
        }

        for (int col = left; col <= right; col++) {
            screen.fillColor(t.color, col, row);
        }
    }
}

int main()
{
    vtkImageData *image = NewImage(1786, 1344);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = 1786*1344;
    for (int i = 0; i < npixels*3; i++) {
        buffer[i] = 0;
    }

    std::vector<Triangle> triangles = GetTriangles();

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1786;
    screen.height = 1344;

    for (int i = 0; i < triangles.size(); i++) {
        triangleID = i;
        if (isGoingDown(triangles[i])) {
            RasterizeGoingDownTriangle(triangles[i], screen);
        }

        if (isGoingUp(triangles[i])) {
            RasterizeGoingUpTriangle(triangles[i], screen);
        }

        if (isArbitrary(triangles[i])) {
            RasterizeArbitraryTriangle(triangles[i], screen);
        }
    }
    WriteImage(image, "allTriangles");
}


