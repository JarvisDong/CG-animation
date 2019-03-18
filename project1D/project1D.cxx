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

int triangleID = -1;
class Triangle
{
public:
    double X[3];
    double Y[3];
    double Z[3];
    double colors[3][3];
};

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}


class Screen
{
public:
    double          width;
    double          height;
    unsigned char   *buffer;
    double          *depthBuffer;
    void fillColor(double color_rc[3], int c, int r, double z) {
        if (r >= height || c < 0 || c >= width || r < 0) { return; }
        int index = (r*width+c)*3;
        if (z > depthBuffer[index]) {
            for (int i = 0; i < 3; i++) {
                buffer[index+i] = ceil_441(color_rc[i]*255);
            }
            depthBuffer[index] = z;
        }

    }
};

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
    rdr->SetFileName("proj1d_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
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
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

double Lerp(double X, double field_A, double field_B, double A, double B) {
    // return the field value of X 
    if (X == A) { return field_A; }
    return field_A+(((X-A)/(B-A))*(field_B-field_A));
}

void RasterizeGoingDownTriangle(Triangle &t, Screen screen) {
    double x_left = 0.0;
    double x_right = 0.0;
    double x_down = t.X[0];
    double y_up_left = 0.0;
    double y_up_right = 0.0;
    double y_down = t.Y[0];
    double z_down = t.Z[0];
    double z_left, z_right;
    double *color_down = t.colors[0];
    double *color_left, *color_right;

    for (int i=1; i<3; i++) {
        if (t.Y[i] < y_down) {
          y_down = t.Y[i];
          x_down = t.X[i];
          z_down = t.Z[i];
          color_down = t.colors[i];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            if (t.Y[i] == t.Y[j]) {
                if (t.X[j] >= t.X[i]) {
                    x_left = t.X[i];
                    y_up_left = t.Y[i];
                    z_left = t.Z[i];
                    color_left = t.colors[i];

                    x_right = t.X[j];
                    y_up_right = t.Y[j];
                    z_right = t.Z[j];
                    color_right = t.colors[j];
                }
                else {
                    x_left = t.X[j];
                    y_up_left = t.Y[j];
                    z_left = t.Z[j];
                    color_left = t.colors[j];

                    x_right = t.X[i];
                    y_up_right = t.Y[i];
                    z_right = t.Z[i];
                    color_right = t.colors[i];
                }
            }
        }
    }
    // printf("working on triangle %d\nRasterizing GoingDownTriangle.\nTriangle:\n", triangleID);
    // printf("Vertex 0: pos = (%f, %f, %f), color = (%f, %f, %f)\n", x_left, y_up_left, z_left, color_left[0], color_left[1], color_left[2]);
    // printf("Vertex 1: pos = (%f, %f, %f), color = (%f, %f, %f)\n", x_down, y_down, z_down, color_down[0], color_down[1], color_down[2]);
    // printf("Vertex 2: pos = (%f, %f, %f), color = (%f, %f, %f)\n", x_right, y_up_right, z_right, color_right[0], color_right[1], color_right[2]);

    /*
    Determine all 3 vertices position for going down triangle
        (x_left, y_up_left, z_left, *color_left)     (x_right, y_up_right, z_right, *color_right)
                                  *  * *  *  *  *  *  *  *  * *  *                                                        
                                    *                           *
                                      *                       *
                                        *                   *
                                          *               *
                                            *           *
                                              *       *
                                                *   *
                                                  * 
                                (x_down, y_down, z_down, *color_down)
    */

    for (int r = ceil_441(y_down); r <= floor_441(y_up_left); r++) {
/*
        double Lerp(double X, double field_A, double field_B, double A, double B) {
            // return the field value of X 
            if (X == A) { return field_A; }
            return field_A+(((X-A)/(B-A))*(field_B-field_A));
        }
*/      
//      Interpolate x(leftEnd) and x(rightEnd) from triangle vertices y values, field = x
        double leftEnd = Lerp(r, x_down, x_left, y_down, y_up_left);
        double rightEnd = Lerp(r, x_down, x_right, y_down, y_up_right);

//      Interpolate z(leftEnd) and z(rightEnd) from triangle vertices y values, field = z

        double Z_leftEnd = Lerp(r, z_left, z_down, y_up_left, y_down);
        double Z_rightEnd = Lerp(r, z_right, z_down, y_up_right, y_down);

//      Calculate Color(leftEnd) and Color(rightEnd) using interpolation from triangle vertices
        double color_leftEnd[3], color_rightEnd[3];
        for (int i = 0; i < 3; i++) {
            color_leftEnd[i] = Lerp(r, color_down[i], color_left[i], y_down, y_up_left);
            color_rightEnd[i] = Lerp(r, color_down[i], color_right[i], y_down, y_up_right);
        }
        // printf("    Rasterizing along row %d with left end %f (Z: %f, RGB = %f/%f/%f) and right end = %f (Z: %f, RGB = %f/%f/%f)\n", r, leftEnd, Z_leftEnd, color_leftEnd[0], color_leftEnd[1], color_leftEnd[2], rightEnd, Z_rightEnd,  color_rightEnd[0],color_rightEnd[1],color_rightEnd[2]);
        
//      Interpolate z(r,c) from z(leftEnd) and z(rightEnd)
        double Z_rc, color_rc[3];
        for (int c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++) {
            Z_rc = Lerp(c, Z_leftEnd, Z_rightEnd, leftEnd, rightEnd);
            for (int i = 0; i < 3; i++) {
                color_rc[i] = Lerp(c, color_leftEnd[i], color_rightEnd[i], leftEnd, rightEnd);
            }
            screen.fillColor(color_rc, c, r, Z_rc);
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
    double z_up = t.Z[0]; 
    double z_left, z_right;
    double *color_up = t.colors[0]; 
    double *color_left, *color_right;

    for (int i=1; i<3; i++) {
        if (t.Y[i] > y_up) {
          y_up = t.Y[i];
          x_up = t.X[i];
          z_up = t.Z[i];
          color_up = t.colors[i];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            if (t.Y[i] == t.Y[j]) {
                if (t.X[j] >= t.X[i]) {
                    x_left = t.X[i];
                    y_down_left = t.Y[i];
                    z_left = t.Z[i];
                    color_left = t.colors[i];

                    x_right = t.X[j];
                    y_down_right = t.Y[j];
                    z_right = t.Z[j];
                    color_right = t.colors[j];
                }
                else {
                    x_left = t.X[j];
                    y_down_left = t.Y[j];
                    z_left = t.Z[j];
                    color_left = t.colors[j];

                    x_right = t.X[i];
                    y_down_right = t.Y[i];
                    z_right = t.Z[i];
                    color_right = t.colors[i];
                }
            }
        }
    }
    // printf("working on triangle %d\nRasterizing GoingUpTriangle.\nTriangle:\n", triangleID);
    // printf("Vertex 0: pos = (%f, %f, %f), color = (%f, %f, %f)\n", x_left, y_down_left, z_left, color_left[0], color_left[1], color_left[2]);
    // printf("Vertex 1: pos = (%f, %f, %f), color = (%f, %f, %f)\n", x_up, y_up, z_up, color_up[0], color_up[1], color_up[2]);
    // printf("Vertex 2: pos = (%f, %f, %f), color = (%f, %f, %f)\n", x_right, y_down_right, z_right, color_right[0], color_right[1], color_right[2]);
    /*     Determine all 3 vertices position for going up triangle

                          (x_up, y_up, z_up, *color_up)
                                        *
                                      *   *
                                    *       *
                                  *           *
                                *               *
                              *                   *
                            * * * * * * * * * * * * *
        (x_left, y_down_left,                       (x_right, y_down_right,
        z_left, *color_left )                       z_right, *color_right )
    */                                              

    for (int r = ceil_441(y_down_right); r <= floor_441(y_up); r++) {
//      Interpolate x(leftEnd) and x(rightEnd) from triangle vertices y values, field = x
        double leftEnd = Lerp(r, x_up, x_left, y_up, y_down_left);
        double rightEnd = Lerp(r, x_up, x_right, y_up, y_down_right);

//      Interpolate z(leftEnd) and z(rightEnd) from triangle vertices y values, field = z
        double Z_leftEnd = Lerp(r, z_up, z_left, y_up, y_down_left);
        double Z_rightEnd = Lerp(r, z_up, z_right, y_up, y_down_right);

//      Calculate Color(leftEnd) and Color(rightEnd) using interpolation from triangle vertices
        double color_leftEnd[3], color_rightEnd[3];
        for (int i = 0; i < 3; i++) {
            color_leftEnd[i] = Lerp(r, color_up[i], color_left[i], y_up, y_down_left);
            color_rightEnd[i] = Lerp(r, color_up[i], color_right[i], y_up, y_down_right);
        }
    
        // printf("    Rasterizing along row %d with left end %f (Z: %f, RGB = %f/%f/%f) and right end = %f (Z: %f, RGB = %f/%f/%f)\n", r, leftEnd, Z_leftEnd, color_leftEnd[0], color_leftEnd[1], color_leftEnd[2], rightEnd, Z_rightEnd,  color_rightEnd[0],color_rightEnd[1],color_rightEnd[2]);


//      Interpolate z(r,c) from z(leftEnd) and z(rightEnd)
        double Z_rc, color_rc[3];
        
        for (int c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++) {
            Z_rc = Lerp(c, Z_leftEnd, Z_rightEnd, leftEnd, rightEnd);
            for (int i = 0; i < 3; i++) {
                color_rc[i] = Lerp(c, color_leftEnd[i], color_rightEnd[i], leftEnd, rightEnd);
            }
            screen.fillColor(color_rc, c, r, Z_rc);
        }
    }
}


void RasterizeArbitraryTriangle(Triangle &t, Screen screen)
{

//  relocate vertices
    double y_down = t.Y[0];
    double y_up = t.Y[0];
    double x_down = t.X[0];
    double x_up = t.X[0];
    double x_mid = 0.0;
    double y_mid = 0.0;
    double z_up = t.Z[0];
    double *color_up = t.colors[0];
    double z_down = t.Z[0];
    double *color_down = t.colors[0];
    double z_mid, *color_mid;


    for (int i = 1; i < 3; i++) {
        if (t.Y[i] > y_up) { 
            y_up = t.Y[i];
            x_up = t.X[i];
            z_up = t.Z[i];
            color_up = t.colors[i];
        }
        if (t.Y[i] < y_down) {
            y_down = t.Y[i];
            x_down = t.X[i];
            z_down = t.Z[i];
            color_down = t.colors[i];
        }
    }

    for (int i = 0; i < 3; i++) {
        if ((t.Y[i] < y_up) && (t.Y[i] > y_down)) {
            y_mid = t.Y[i];
            x_mid = t.X[i];
            z_mid = t.Z[i];
            color_mid = t.colors[i];
        }
    }

    // double x_mid_alt = Lerp(y_mid, x_up, x_down, y_up, y_down);
    double slope_up_down = (y_up - y_down)/(x_up-x_down);
    double x_mid_alt = (y_mid-y_up)/slope_up_down + x_up;
    double y_mid_alt = y_mid;
    double z_mid_alt = Lerp(y_mid, z_up, z_down, y_up, y_down);
    double color_mid_alt[3];
    for (int i = 0; i < 3; i++) {
        color_mid_alt[i] = Lerp(y_mid_alt, color_up[i], color_down[i], y_up, y_down);
    }

/*     Determine all 3 vertices position for arbitrary triangle
                                        (x_up, y_up, z_up, *color_up)       
                                                    *
                                                  * *
                                                *   *
                                              *     *
                                            *       *
       (x_mid, y_mid, z_mid, *color_mid)  *---------* (x_mid_alt, y_mid, z_mid_alt, *color_mid_alt)
                                            *       *
                                              *     *
                                                *   *
                                                  * *
                                                    *     
                                        (x_down, y_down, z_down, *color_down)  
*/
// Going up
    // printf("working on triangle %d\nRasterizing GoingUpTriangle.\nTriangle:\n", triangleID);
    // printf("Vertex 0: pos = (%f, %f, %f), color = (%f, %f, %f)\n", x_mid_alt, y_mid, z_mid_alt, color_mid_alt[0], color_mid_alt[1], color_mid_alt[2]);
    // printf("Vertex 1: pos = (%f, %f, %f), color = (%f, %f, %f)\n", x_mid, y_mid, z_mid, color_mid[0], color_mid[1], color_mid[2]);
    // printf("Vertex 2: pos = (%f, %f, %f), color = (%f, %f, %f)\n", x_up, y_up, z_up, color_up[0], color_up[1], color_up[2]);

    Triangle top, bottom;
    top.X[0] = x_up;
    top.X[1] = x_mid;
    top.X[2] = x_mid_alt;
    top.Y[0] = y_up;
    top.Y[1] = y_mid;
    top.Y[2] = y_mid_alt;
    top.Z[0] = z_up;
    top.Z[1] = z_mid;
    top.Z[2] = z_mid_alt;

    bottom.X[0] = x_down;
    bottom.X[1] = x_mid;
    bottom.X[2] = x_mid_alt;
    bottom.Y[0] = y_down;
    bottom.Y[1] = y_mid;
    bottom.Y[2] = y_mid_alt;
    bottom.Z[0] = z_down;
    bottom.Z[1] = z_mid;
    bottom.Z[2] = z_mid_alt;

    for (int i = 0; i < 3; i++) {
        top.colors[0][i] = color_up[i];
        top.colors[1][i] = color_mid[i];
        top.colors[2][i] = color_mid_alt[i];
        bottom.colors[0][i] = color_down[i];
        bottom.colors[1][i] = color_mid[i];
        bottom.colors[2][i] = color_mid_alt[i];
    }
    RasterizeGoingUpTriangle(top, screen);
    RasterizeGoingDownTriangle(bottom, screen);
}

int main()
{
    vtkImageData *image = NewImage(1000, 1000);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = 1000*1000;
    double *depthBuffer = new double[npixels*3];
    for (int i = 0; i < npixels*3; i++) {
        buffer[i] = 0;
        depthBuffer[i] = -1;
    }
    
    std::vector<Triangle> triangles = GetTriangles();

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;
    screen.depthBuffer = depthBuffer;
    int size = triangles.size();
    for (int i = 0; i < size; i++) {
        triangleID = i;
        Triangle T = triangles[i];
        if (T.Y[0] == T.Y[1]){ 
            if (T.Y[2] < T.Y[0]){ RasterizeGoingDownTriangle(T, screen); }
            else{ RasterizeGoingUpTriangle(T, screen); }
        }
        else if (T.Y[0] == T.Y[2]) {
            if (T.Y[1] < T.Y[0]){ RasterizeGoingDownTriangle(T, screen); }
            else{ RasterizeGoingUpTriangle(T, screen); }
        }
        else if (T.Y[1] == T.Y[2]){
            if (T.Y[0] < T.Y[1]){ RasterizeGoingDownTriangle(T, screen); }
            else{ RasterizeGoingUpTriangle(T, screen); }
        }
        else { RasterizeArbitraryTriangle(T, screen); }
    }
    WriteImage(image, "allTriangles");
}
