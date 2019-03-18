#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>

class Triangle
{
public:
    double X[3];
    double Y[3];
    double Z[3];
    double colors[3][3];
};
int triangleID = 0;
double ceiled(double f)
{
    return ceil(f-0.00001);
}

double floored(double f)
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
    void fillColor(double color_rc[3], int c, int r) {
        if (r >= height || c < 0 || c >= width || r < 0) { return; }
        int index = (r*width+c)*3;
        for (int i = 0; i < 3; i++) {
            buffer[index+i] = ceiled(color_rc[i]*255);
        }
    }
};


class Matrix
{
  public:
    double          A[4][4];        // A[i][j] means row i, column j
    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

class Camera
{
  public:
    double   near, far;                  // Distance between camera with near/far plane
    double   angle;                      // in radian, Cmath work with radian
    double   position[3];                // Position of the camera
    double   focus[3];                   // focu point
    double   up[3];                      /* direction from the base of "nose" to "forehead", eyes as camera
                                            (determine which triangle should on bottom and which should be on top)
                                            the direction is always perpendicular to the line camera position to focu point. 
                                         */
    // double   frame[4][3];                                    
    Matrix   CameraTransform(void);         // World space  -> Camera space return a Matrix CT
    Matrix   ViewTransform(void);           // Camera space -> Image space  return a Matrix VT
    Matrix   DeviceTransform(void);         // Image space -> Device space retrun a Matrix DT
};


std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
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
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
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
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

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

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    /*
    Not use for project1E, but will be useful for porject1F
    */
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
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

Matrix
Camera::ViewTransform() {
/*
    input parameters: (α,n,f)
    image space cube -1 <= u,v,w <= 1
    [cos(α/2)   0           0             0]
    [0          cos(α/2)    0             0]
    [0          0           (f+n)/(f-n)  -1]
    [0          0           2fn/(f-n)     0]
*/
    Matrix CameraMatrix;
    CameraMatrix.A[0][0] = 1/(tan(angle/2));
    CameraMatrix.A[0][1] = 0;
    CameraMatrix.A[0][2] = 0; 
    CameraMatrix.A[0][3] = 0;

    CameraMatrix.A[1][0] = 0;
    CameraMatrix.A[1][1] = 1/(tan(angle/2));
    CameraMatrix.A[1][2] = 0;
    CameraMatrix.A[1][3] = 0;

    CameraMatrix.A[2][0] = 0; 
    CameraMatrix.A[2][1] = 0;
    CameraMatrix.A[2][2] = (far+near)/(far-near);
    CameraMatrix.A[2][3] = -1;

    CameraMatrix.A[3][0] = 0;
    CameraMatrix.A[3][1] = 0;
    CameraMatrix.A[3][2] = (2*far*near)/(far-near);
    CameraMatrix.A[3][3] = 0;

    return CameraMatrix;
}

Matrix
Camera::DeviceTransform() {
    Matrix imageMatrix;
/*
                            (n/2  0    0   0)
    (x  y  z  1)    *       (0    m/2  0   0)
                            (0    0    1   0)
                            (n/2  m/2  0   1)
*/
    imageMatrix.A[0][0] = 1000/2;  imageMatrix.A[0][1] = 0;      imageMatrix.A[0][2] = 0;  imageMatrix.A[0][3] = 0;
    imageMatrix.A[1][0] = 0;       imageMatrix.A[1][1] = 1000/2; imageMatrix.A[1][2] = 0;  imageMatrix.A[1][3] = 0;
    imageMatrix.A[2][0] = 0;       imageMatrix.A[2][1] = 0;      imageMatrix.A[2][2] = 1;  imageMatrix.A[2][3] = 0;
    imageMatrix.A[3][0] = 1000/2;  imageMatrix.A[3][1] = 1000/2; imageMatrix.A[3][2] = 0;  imageMatrix.A[3][3] = 1;

    return imageMatrix;    
}

Matrix 
Camera::CameraTransform() {
    double O[3], w[3], v[3], u[3];
    /*
        Construct a camera frame
        * must choose (u, v, w, O)
        * O = camera position 
        * u = Up x (O-focus) 
        * v = (O-focus) x u 
        * w = O-focus
        * note: u, v need cross product
    */
    for (int i = 0; i < 3; i++) {
        O[i] = position[i];
        w[i] = O[i] - focus[i];
    }
    /* Cross product for u, v
    AxB = (A.y*B.z - A.z*B.y,
           A.z*B.x - A.x*B.z,
           A.x*B.y - A.y*B.x)
    */
    u[0] = up[1]*w[2] - up[2]*w[1];
    u[1] = up[2]*w[0] - up[0]*w[2];
    u[2] = up[0]*w[1] - up[1]*w[0];

    v[0] = w[1]*u[2] - w[2]*u[1];
    v[1] = w[2]*u[0] - w[0]*u[2];
    v[2] = w[0]*u[1] - w[1]*u[0];

    double u_norm = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
    double v_norm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    double w_norm = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);

    for (int i = 0; i < 3; i++) {
        u[i] = u[i]/u_norm;
        v[i] = v[i]/v_norm;
        w[i] = w[i]/w_norm;
    }
    
    /*
        Construct a matrix to transform points
        from Cartesian Frame to Camera Frame

        [u.x    v.x     w.x     O]

        [u.y    v.y     w.y     O]

        [u.z    v.z     w.z     O]

        [u*t    v*t     w*t     1]

        where t = (0, 0, 0) - O
    */
    Matrix deviceMatrix;
    double t[3] = {0-O[0], 0-O[1], 0-O[2]};

    for (int i = 0; i < 3; i++) {
        deviceMatrix.A[i][0] = u[i];
        deviceMatrix.A[i][1] = v[i];
        deviceMatrix.A[i][2] = w[i];
        deviceMatrix.A[i][3] = 0;   
    }
    deviceMatrix.A[3][0] = u[0]*t[0] + u[1]*t[1] + u[2]*t[2];
    deviceMatrix.A[3][1] = v[0]*t[0] + v[1]*t[1] + v[2]*t[2];
    deviceMatrix.A[3][2] = w[0]*t[0] + w[1]*t[1] + w[2]*t[2];
    deviceMatrix.A[3][3] = 1;

    return deviceMatrix;
}

double Lerp(double X, double field_A, double field_B, double A, double B) {
    // return the field value of X
    if (X == A) { return field_A; }
    // if (floored(fabs(X - A)) == 0) { return field_A; }
    // if (floored(fabs(B - A)) == 0) { return field_A; }
    return field_A+(((X-A)/(B-A))*(field_B-field_A));
}

double GetIntersectionValueX(int scanLine, double x1, double x2, double y1, double y2)
{
    double doublescanline = (double)scanLine;
    if (floored(fabs(x1 - x2)) == 0) { return x1; }
    else {
        double slope = (y1-y2)/(x1-x2);
        if (floored(fabs(y1 - y2)) == 0) { return x1; }
        return (doublescanline-y1)/slope + x1;
    }
}

void RasterizeGoingDownTriangle(Triangle &t, Screen screen) {
    double x_left = t.X[0];
    double x_right = t.X[0];
    double x_down = t.X[0];
    double y_up_left = t.Y[0];
    double y_up_right = t.Y[0];
    double y_down = t.Y[0];
    double z_down = t.Z[0];
    double z_left = t.Z[0];
    double z_right = t.Z[0];
    double *color_down = t.colors[0];
    double *color_left = t.colors[0];
    double *color_right = t.colors[0];

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



    /*
    Determine all 3 vertices position for going down triangle
                vertex    A                                      vertex     B
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
                                            vertex  C
    */

    for (int r = ceiled(y_down); r <= floored(y_up_left); r++) {
 
//      Interpolate x(leftEnd) and x(rightEnd) from triangle vertices y values, field = x
        double leftEnd = Lerp(r, x_left, x_down, y_up_left, y_down);
        double rightEnd = Lerp(r, x_right, x_down, y_up_right, y_down);
        // double leftEnd = GetIntersectionValueX(r, x_left, x_down, y_up_left, y_down);
        // double rightEnd = GetIntersectionValueX(r, x_right, x_down, y_up_right, y_down);

//      Interpolate z(leftEnd) and z(rightEnd) from triangle vertices y values, field = z

        double Z_leftEnd = Lerp(r, z_down, z_left, y_up_left, y_down);
        double Z_rightEnd = Lerp(r, z_down, z_right, y_up_right, y_down);

//      Calculate Color(leftEnd) and Color(rightEnd) using interpolation from triangle vertices
        double color_leftEnd[3], color_rightEnd[3];
        for (int i = 0; i < 3; i++) {
            color_leftEnd[i] = Lerp(r, color_down[i], color_left[i], y_down, y_up_left);
            color_rightEnd[i] = Lerp(r, color_down[i], color_right[i], y_down, y_up_right);
        }
        // printf("    Rasterizing along row %d with left end %f (Z: %f, RGB = %f/%f/%f) and right end = %f (Z: %f, RGB = %f/%f/%f)\n", r, leftEnd, Z_leftEnd, color_leftEnd[0], color_leftEnd[1], color_leftEnd[2], rightEnd, Z_rightEnd,  color_rightEnd[0],color_rightEnd[1],color_rightEnd[2]);
        
//      Interpolate z(r,c) from z(leftEnd) and z(rightEnd)
        double z, color_rc[3];
        for (int c = ceiled(leftEnd); c <= floored(rightEnd); c++) {
            z = Lerp(c, Z_leftEnd, Z_rightEnd, leftEnd, rightEnd);
            for (int i = 0; i < 3; i++) {
                color_rc[i] = Lerp(c, color_leftEnd[i], color_rightEnd[i], leftEnd, rightEnd);
            }
            int index = (r*1000+c)*3;
            if (z > screen.depthBuffer[index]) {
                screen.fillColor(color_rc, c, r);
                // printf("            Got fragment r = %d, c = %d, z = %f, color = %f/%f/%f\n", r, c, Z_rc, color_rc[0], color_rc[1], color_rc[2]);
                screen.depthBuffer[index] = z;
            }            
        }
    }
}

void RasterizeGoingUpTriangle(Triangle &t, Screen screen) {
    double x_up = t.X[0];
    double x_left = t.X[0];
    double x_right = t.X[0];
    double y_up = t.Y[0];
    double y_down_left = t.Y[0];
    double y_down_right = t.Y[0];
    double z_up = t.Z[0]; 
    double z_left = t.Z[0];
    double z_right = t.Z[0];
    double *color_up = t.colors[0]; 
    double *color_left = t.colors[0];
    double *color_right = t.colors[0];

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

    /*     Determine all 3 vertices position for going up triangle
                                    vertex  A
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
            vertex  B                                       vertex C
    */                                              

    for (int r = ceiled(y_down_right); r <= floored(y_up); r++) {
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
        double z, color_rc[3];
        
        for (int c = ceiled(leftEnd); c <= floored(rightEnd); c++) {
            z = Lerp(c, Z_leftEnd, Z_rightEnd, leftEnd, rightEnd);
            for (int i = 0; i < 3; i++) {
                color_rc[i] = Lerp(c, color_leftEnd[i], color_rightEnd[i], leftEnd, rightEnd);
            }
            int index = (r*1000+c)*3;
            if (z > screen.depthBuffer[index]) {
                screen.fillColor(color_rc, c, r);
                // printf("            Got fragment r = %d, c = %d, z = %f, color = %f/%f/%f\n", r, c, Z_rc, color_rc[0], color_rc[1], color_rc[2]);
                screen.depthBuffer[index] = z;
            }            
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
    double z_mid = t.Z[0];
    double *color_mid = t.colors[0];


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
    double slope_up_down = (y_up - y_down)/(x_up-x_down);
    double x_mid_alt = (y_mid-y_up)/slope_up_down + x_up;

    double y_mid_alt = y_mid;
    double z_mid_alt = Lerp(y_mid, z_up, z_down, y_up, y_down);

    double color_mid_alt[3];
    for (int i = 0; i < 3; i++) {
        color_mid_alt[i] = Lerp(x_mid_alt, color_up[i], color_down[i], x_up, x_down);
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

int f = 0;
int main() {
    vtkImageData *image = NewImage(1000, 1000);
    std::vector<Triangle> triangles = GetTriangles();

    Screen screen;

    for (int n = 0; n < 1; n++) {
        f = 250*n;  // f = 0, 250, 500, 750
        char temp[32];
        sprintf(temp, "frame%d", f);

        /* 
            **** initialize screen ****
        need to initialize the color and z-buffers for each rendering.
        */
        unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
        int npixels = 1000*1000;
        double *depthBuffer = new double[npixels*3];
        for (int i = 0; i < npixels*3; i++) {
            buffer[i] = 0;
            depthBuffer[i] = -1;
        }

        screen.width = 1000;
        screen.height = 1000;
        screen.buffer = buffer;
        screen.depthBuffer = depthBuffer;

        Camera c = GetCamera(f, 1000); // called 4 times

        // compose 3 transforms into one. Finally we want CTxVTxDT
        Matrix CT,VT,DT,VTxDT,M;
        CT = c.CameraTransform();
        VT = c.ViewTransform();
        DT = c.DeviceTransform();
        VTxDT = VTxDT.ComposeMatrices(VT, DT);
        M = M.ComposeMatrices(CT, VTxDT);
        /*
         ******** Transform triangles to Device space ********
            involves setting up and applying matrices
            if vector<Triangle> t get modified, remember
            to undo it later
        */

        int size = triangles.size();
        for (int i = 0; i < size; i++) {
            triangleID = i;
            Triangle t = triangles[i];
            double * Ain = new double[4];
            double * Bin = new double[4];
            double * Cin = new double[4];
            double * Aout = new double[4];
            double * Bout = new double[4];
            double * Cout = new double[4];

            Ain[0] = t.X[0];  Ain[1] = t.Y[0];  Ain[2] = t.Z[0];  Ain[3] = double(1);
            Bin[0] = t.X[1];  Bin[1] = t.Y[1];  Bin[2] = t.Z[1];  Bin[3] = double(1);
            Cin[0] = t.X[2];  Cin[1] = t.Y[2];  Cin[2] = t.Z[2];  Cin[3] = double(1);

            M.TransformPoint(Ain, Aout);
            M.TransformPoint(Bin, Bout);
            M.TransformPoint(Cin, Cout);

            Triangle T;
            T.X[0] = Aout[0]/Aout[3];
            T.X[1] = Bout[0]/Bout[3];
            T.X[2] = Cout[0]/Cout[3];
            T.Y[0] = Aout[1]/Aout[3];
            T.Y[1] = Bout[1]/Bout[3];
            T.Y[2] = Cout[1]/Cout[3];
            T.Z[0] = Aout[2]/Aout[3];
            T.Z[1] = Bout[2]/Bout[3];
            T.Z[2] = Cout[2]/Cout[3];
            for (int j = 0; j < 3; j++) {
                T.colors[0][j] = t.colors[0][j];
                T.colors[1][j] = t.colors[1][j];
                T.colors[2][j] = t.colors[2][j];
            }

        // Rasterize
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
        // Save images
        if (n == 0){
            WriteImage(image, "frame000");
        }
        else{
            WriteImage(image, temp);
        }
    }
}