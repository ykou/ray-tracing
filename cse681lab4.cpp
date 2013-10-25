//
//  cse681lab4.cpp
//  cse681lab4
//
//  Created by Kou Yuxiang on 12/2/11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector.h>

#include <Inventor/SbLinear.h>
#include <Inventor/actions/SoGetMatrixAction.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoLight.h>
#include <Inventor/nodes/SoPointLight.h>

#include "OSUInventor.h"

#define MAXRECURSION 3
#define AIRINDEX 1.0
#define SPHEREINDEX 1.5
#define PI 3.14159

using namespace std;

int shadowOn = 0;
int reflectionOn = 0;
int superSampleOn = 0;
int depthOfFieldOn = 0;
float randomNo();

class cameraClass {
public:
    SoCamera * camera;
    SbVec3f position, direction, viewUp, rotationAxis, 
    n, u, v;
    SbRotation orientation;
    float rotationAngle, aspectRatio, viewAngle, 
    W, H ;
};

class pointLightClass {
public:
    SbVec3f location, color;
    float intensity;
};

class sphereClass {
public: 
    int type;   // 0 is sphere, 1 is cube, 2 is TODO
    SbVec3f center, ambientColor, diffuseColor, specularColor;
    float radius, width, height, depth;
    float shininess, transparency;
    SbMatrix F;
};

class rayTracerClass {
public:
    sphereClass *mySphere;
    cameraClass myCamera;
    pointLightClass *myPointLight;
    int resolutionX, resolutionY, objectNumber, lightNumber;
    float *image;
    
    void readOpenInventorScene(int argc, char *argv, 
                               char *resX, char *resY);
    void traceRays(int imageHeight, int imageWidth);
    bool intersect(SbVec3f rayStart, SbVec3f rayEnd);
    SbVec3f shade(SbVec3f viewOri, SbVec3f viewDir, int recursionDep);
};

void usageError();
SbVec3f scaleColor(SbVec3f sourceColor);
SbVec3f backgroundColor(SbVec3f sourceColor);
SbVec3f white(1, 1, 1);

SbVec3f shade(float* rayOrigin, float* rayDirection, SbVec3f color);

int main(int argc, char **argv) {
    time_t seconds; 
    time(&seconds);  // this function returns the current time to the variable seconds
    srand((unsigned int)seconds); 
    
    rayTracerClass myRayTracer;
    myRayTracer.readOpenInventorScene(argc, argv[1], argv[2], argv[3]);
    superSampleOn = atoi(argv[4]);
    shadowOn = atoi(argv[5]);
    reflectionOn = atoi(argv[6]);
    depthOfFieldOn = atoi(argv[7]);
    cout << "shadonw On is: " << shadowOn << endl;
    cout << "reflection On is: " << reflectionOn << endl;
    myRayTracer.traceRays(myRayTracer.resolutionX, myRayTracer.resolutionY);
    cout << "random: " << randomNo() << endl;
    cout << "random: " << randomNo() << endl;
    cout << "random: " << randomNo() << endl;
}

void rayTracerClass::readOpenInventorScene(int argc, char *argv, char *resX, char *resY) {
    if (argc != 8)
        usageError();
    
    SoDB::init();
    OSUInventorScene *scene = new OSUInventorScene(argv);
    objectNumber = scene->Objects.getLength();
    mySphere = new sphereClass[objectNumber];
    
    cout << "Read objects from file " << argv << "." << endl;
    cout << "Number of objects = " 
    << objectNumber << "." << endl;
    
    // list objects
    for (int i = 0; i < objectNumber; i++) {
        OSUObjectData * obj = (OSUObjectData *)scene->Objects[i];
        
        if (!obj->Check()) {
            cerr << "Error detected in OSUObjectData for object " << i << "." 
            << endl;
            exit(20);
        };
        
        SoType shape_type = obj->shape->getTypeId();
        
        cout << "Object " << i << " is a "
        << shape_type.getName().getString() << "." << endl;
        
        if (shape_type == SoSphere::getClassTypeId()) {
            SoSphere * sphere = (SoSphere *) obj->shape;
            mySphere[i].type = 0;
            mySphere[i].radius = sphere->radius.getValue();
            cout << "  Sphere radius = " << mySphere[i].radius << "."
            << endl;
        };
        
        if (shape_type == SoCube::getClassTypeId()) {
            SoCube * cube = (SoCube *) obj->shape;
            mySphere[i].type = 1;
            mySphere[i].width = cube->width.getValue();
            mySphere[i].height = cube->height.getValue();
            mySphere[i].depth = cube->depth.getValue();
            cout << "  Cube width, height, depth = (" 
            << mySphere[i].width << ","
            << mySphere[i].height << ","
            << mySphere[i].depth << ")." << endl;
        };
        
        // object transformation
        SoTransform * transformation = obj->transformation;
        SbVec3f scale_vector = transformation->scaleFactor.getValue();
        scale_vector[0] = 0.5;
        cout << "  scale_vector = ("
        << scale_vector[0] << ","
        << scale_vector[1] << ","
        << scale_vector[2] << ")" << endl;
        SbRotation rotation = transformation->rotation.getValue();
        SbVec3f rotationAxis;
        float rotationAngle;
        rotation.getValue(rotationAxis, rotationAngle);
        SbVec3f translation_vector = transformation->translation.getValue();
        
        mySphere[i].center[0] = 0;
        mySphere[i].center[1] = 0;
        mySphere[i].center[2] = 0;
        SbMatrix T, S, R, F; 
        T.setTranslate(translation_vector); 
        R.setRotate(rotation); 
        S.setScale(scale_vector);
        F = S * R * T;
        mySphere[i].F = F;
        //SbMatrix TEIQWEOI(1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
        //F = F.multRight(TEIQWEOI);
        //TEIQWEOI.multVecMatrix(mySphere[i].center, mySphere[i].center);
        F.multVecMatrix(mySphere[i].center, mySphere[i].center);
        if (mySphere[i].type == 0) {
            //mySphere[i].radius = mySphere[i].radius * scale_vector[0];
            cout << "  after translation Sphere radius:" << mySphere[i].radius << endl;
        }
        if (mySphere[i].type == 1) {
            cout << "  after translation Cube width, height, depth = (" 
            << mySphere[i].width << ","
            << mySphere[i].height << ","
            << mySphere[i].depth << ")." << endl;
        }
        
        cout << "  center is ("
        << mySphere[i].center[0] << ","
        << mySphere[i].center[1] << ","
        << mySphere[i].center[2] << ")." << endl;
        
        cout << "  Scale by (" 
        << scale_vector[0] << ","
        << scale_vector[1] << ","
        << scale_vector[2] << ")." << endl;
        cout << "  Rotate around axis (" 
        << rotationAxis[0] << ","
        << rotationAxis[1] << ","
        << rotationAxis[2] << ")" 
        << " by " << rotationAngle << " radians." << endl;
        cout << "  Translate by (" 
        << translation_vector[0] << ","
        << translation_vector[1] << ","
        << translation_vector[2] << ")." << endl;
        
        // object material (color)
        SoMaterial * material = obj->material;
        mySphere[i].transparency = material->transparency[0];
        cout << "  Transparency is: " << mySphere[i].transparency << endl;
        mySphere[i].shininess = material->shininess[0];
        cout << "  Shininess is: " << mySphere[i].shininess << endl;
        // material->diffuseColor[0] is the first entry in an array of colors
        mySphere[i].ambientColor[0] = material->ambientColor[0][0];
        mySphere[i].ambientColor[1] = material->ambientColor[0][1];
        mySphere[i].ambientColor[2] = material->ambientColor[0][2];
        cout << "  Ambient color (red, green, blue) = ("
        << mySphere[i].ambientColor[0] << ","
        << mySphere[i].ambientColor[1] << ","
        << mySphere[i].ambientColor[2] << ")." << endl;
        mySphere[i].diffuseColor[0] = material->diffuseColor[0][0];
        mySphere[i].diffuseColor[1] = material->diffuseColor[0][1];
        mySphere[i].diffuseColor[2] = material->diffuseColor[0][2];
        cout << "  Diffuse color (red, green, blue) = ("
        << mySphere[i].diffuseColor[0] << ","
        << mySphere[i].diffuseColor[1] << ","
        << mySphere[i].diffuseColor[2] << ")." << endl;
        mySphere[i].specularColor[0] = material->specularColor[0][0];
        mySphere[i].specularColor[1] = material->specularColor[0][1];
        mySphere[i].specularColor[2] = material->specularColor[0][2];
        cout << "  Specular color (red, green, blue) = ("
        << mySphere[i].specularColor[0] << ","
        << mySphere[i].specularColor[1] << ","
        << mySphere[i].specularColor[2] << ")." << endl;
    };
    
    // list lights
    cout << endl;
    lightNumber = scene->Lights.getLength();
    myPointLight = new pointLightClass[lightNumber];
    for (int j = 0; j < scene->Lights.getLength(); j++) {
        SoLight * light = (SoLight *) scene->Lights[j];
        SoType light_type = light->getTypeId();
        
        if (light->on.getValue() == true)
            cout << "Light " << j << " is on." << endl;
        else
            cout << "Light " << j << " is off." << endl;
        
        if (light_type == SoPointLight::getClassTypeId()) {
            SoPointLight * point_light = (SoPointLight *) light;
            //SbVec3f location = point_light->location.getValue();
            myPointLight[j].location = point_light->location.getValue();
            myPointLight[j].intensity = point_light->intensity.getValue();
            myPointLight[j].color = point_light->color.getValue();
            cout << "  Point light location = (" 
            << myPointLight[j].location[0] << ","
            << myPointLight[j].location[1] << ","
            << myPointLight[j].location[2] << ")." << endl;
            cout << "  Point light intensity = (" 
            << myPointLight[j].intensity << ")." << endl;
            cout << "  Point light color = (" 
            << myPointLight[j].color[0] << ","
            << myPointLight[j].color[1] << ","
            << myPointLight[j].color[2] << ")." << endl;
        };
    };
    
    // camera
    if (scene->Camera == NULL) {
        cout << endl;
        cout << "No camera found." << endl;
    }
    else {
        SoCamera * camera = scene->Camera;
        //SbVec3f camera_position = camera->position.getValue();
        myCamera.position = camera->position.getValue();
        //SbRotation camera_orientation = camera->orientation.getValue();
        myCamera.orientation = camera->orientation.getValue();
        //SbVec3f camera_rotationAxis;
        //float camera_rotationAngle;
        //camera_orientation.getValue(camera_rotationAxis, camera_rotationAngle);
        myCamera.orientation.getValue(myCamera.rotationAxis, 
                                      myCamera.rotationAngle);
        //float camera_aspectRatio = camera->aspectRatio.getValue();
        myCamera.aspectRatio = camera->aspectRatio.getValue();
        SoType camera_type = camera->getTypeId();
        
        // calculate camera direction and camera up direction
        //SbVec3f camera_direction, camera_up;
        //camera_orientation.multVec(SbVec3f(0, 0, -1), camera_direction);
        myCamera.orientation.multVec(SbVec3f(0, 0, -1), myCamera.direction);
        //camera_orientation.multVec(SbVec3f(0, 1, 0), camera_up);
        myCamera.orientation.multVec(SbVec3f(0, 1, 0), myCamera.viewUp);
        
        cout << endl;
        cout << "Camera position = ("
        << myCamera.position[0] << ","
        << myCamera.position[1] << ","
        << myCamera.position[2] << ")." << endl;
        cout << "Camera rotation axis = (" 
        << myCamera.rotationAxis[0] << ","
        << myCamera.rotationAxis[1] << ","
        << myCamera.rotationAxis[2] << ")." << endl;
        cout << "Camera rotation angle = " << myCamera.rotationAngle << " radians."
        << endl;
        cout << "Camera direction = (" 
        << myCamera.direction[0] << ","
        << myCamera.direction[1] << ","
        << myCamera.direction[2] << ")." << endl;
        cout << "Camera up direction = ("
        << myCamera.viewUp[0] << ","
        << myCamera.viewUp[1] << ","
        << myCamera.viewUp[2] << ")." << endl;
        cout << "Camera width/height aspect ratio = "
        << myCamera.aspectRatio << "." << endl;
        
        
        if (camera_type == SoPerspectiveCamera::getClassTypeId()) {
            // camera is a perspective camera
            SoPerspectiveCamera * perspective_camera = 
            (SoPerspectiveCamera *) camera;
            //float camera_height_angle = perspective_camera->heightAngle.getValue();
            myCamera.viewAngle = perspective_camera->heightAngle.getValue();
            cout << "Perspective camera height angle = "
            << myCamera.viewAngle << " radians." << endl;
        };
    };
    
    resolutionX = atoi(resX);
    resolutionY = atoi(resY);
    image = new float[resolutionX * resolutionY * 3];
    cout << "Resolution is: " 
    << resolutionX 
    << " x " 
    << resolutionY
    << endl;
    //return(1);
}

void rayTracerClass::traceRays(int imageWidth, int imageHeight) {
    //calculate the ray origin and ray direction;
    SbVec3f t1 = myCamera.direction;
    t1.negate();    //*caution!* here I inversed the camera
    SbVec3f N = t1;
    SbVec3f n = N/N.length();
    SbVec3f U = myCamera.viewUp.cross(n);
    SbVec3f u = U/U.length();
    SbVec3f v = n.cross(u);
    myCamera.n = n;
    myCamera.u = u;
    myCamera.v = v;
    
    float d = 15;     //assumption by myself
    myCamera.H = tan(myCamera.viewAngle / 2) * d * 2;
    cout << " myCamera.H:" << myCamera.H << endl;
    myCamera.W = myCamera.H * myCamera.aspectRatio;
    cout << " myCamera.W:" << myCamera.W << endl;
    SbVec3f C = myCamera.position - myCamera.n * d;
    SbVec3f L = C - myCamera.u * (myCamera.W/2) 
    + myCamera.v * (myCamera.H/2);    // top-left pixel on the image
    SbVec3f *s = new SbVec3f[imageWidth * imageHeight];
    
    for (int i = 0; i < imageHeight; i++)
        for (int j = 0; j < imageWidth; j++) {
            
            SbVec3f thisPixelColor(0, 0, 0);
            SbVec3f thisPixelColorForSuperSmaple(0, 0, 0);
            for (int k = 0; k < 15 * superSampleOn + 1 ; k ++) {
                // Pixel Supersampling
                s[i * imageWidth + j] = L + myCamera.u * (j + superSampleOn * randomNo()) * myCamera.W/imageWidth - myCamera.v * (i + superSampleOn * randomNo()) * myCamera.H/imageHeight;
                
                SbVec3f thisPixelColorForDepth(0, 0, 0);
                for (int l = 0; l < 50 * depthOfFieldOn + 1; l ++) {
                    //float tempviewo = rand()/(float)(RAND_MAX);
                    float tempviewo = depthOfFieldOn * (randomNo() - 0.5);
                    SbVec3f viewO = myCamera.position + tempviewo * myCamera.u + tempviewo * myCamera.v;
                    SbVec3f viewD = s[i * imageWidth + j] - viewO;
                    
                    thisPixelColorForDepth += shade(viewO, viewD, 0);
                }
                thisPixelColorForDepth = thisPixelColorForDepth / (50 * (float)depthOfFieldOn + 1);
                thisPixelColorForSuperSmaple += thisPixelColorForDepth;
                /*else {
                 
                 SbVec3f viewO = myCamera.position;
                 SbVec3f viewD = s[i * imageWidth + j] - myCamera.position;
                 
                 thisPixelColor += shade(viewO, viewD, 0);
                 }*/
            }
            thisPixelColor = thisPixelColorForSuperSmaple / (15.0 * (float)superSampleOn + 1.0);
            /*else {
             s[i * imageWidth + j] = L + myCamera.u * j * myCamera.W/imageWidth - myCamera.v * i * myCamera.H/imageHeight;
             
             SbVec3f viewO = myCamera.position;
             SbVec3f viewD = s[i * imageWidth + j] - myCamera.position;
             
             thisPixelColor = shade(viewO, viewD, 0);
             }*/
            thisPixelColor = scaleColor(thisPixelColor);
            
            image[(i * imageWidth + j) * 3] = thisPixelColor[0];
            image[(i * imageWidth + j) * 3 + 1] = thisPixelColor[1];
            image[(i * imageWidth + j) * 3 + 2] = thisPixelColor[2];
        }
    
    
    ofstream myfile;
    myfile.open("testppm.ppm",ios::trunc);
    myfile << "P3" << endl;
    myfile << imageWidth << " " << imageHeight << endl;
    myfile << "255" << endl;
    for(int i = 0; i < imageHeight; i++)
        for (int j=0; j < imageWidth; j++) {
            myfile
            << (int)(image[(i * imageWidth + j) * 3] * 255) << " "
            << (int)(image[(i * imageWidth + j) * 3 + 1] * 255) << " "
            << (int)(image[(i * imageWidth + j) * 3 + 2] * 255) << " "
            << endl;
        }
    
    SbVec3f testa(0,0,2);
    cout << "testa:" << testa[0] << testa[1] << testa[2] << endl;
    cout << "testa length:" << testa.length() << endl;
    cout << "testa sqrLength:" << testa.sqrLength() << endl;
    testa.normalize();
    cout << "testa normalized:" << testa[0] << testa[1] << testa[2] << endl;
    testa = testa * 3;
    cout << "testa * 3:" << testa[0] << testa[1] << testa[2] << endl;
    cout << "sqrt(4)=" << sqrt(4) << endl;
    
    cout << "W and H:" 
    << myCamera.W << " " 
    << myCamera.H << endl;
    cout << "N:" << N[0] << N[1] << N[2] <<endl;
    cout << "n:" << n[0] << n[1] << n[2] <<endl;
    cout << "U:" << U[0] << U[1] << U[2] <<endl;
    cout << "u:" << u[0] << u[1] << u[2] <<endl;
    cout << "v:" << v[0] << v[1] << v[2] <<endl;
    
}

void usageError() {
    cerr << "Wrong command!" 
    << endl 
    << "Usage: ./rt {filename} {imageWidth} {imageHeight} {supersampleOn} {shadowOn} {reflectionOn} {depthoffieldOn}" << endl;
    exit(10);
}

SbVec3f rayTracerClass::shade(SbVec3f viewOri, SbVec3f viewDir, int recursionDep) {
    float *t = new float[objectNumber];
    SbVec3f finalColor(0, 0, 0);
    SbVec3f viewOrigin = viewOri;
    SbVec3f viewDirection = viewDir;
    int recursionDepth = recursionDep;
    
    for (int k = 0; k < objectNumber; k++) {
        float tmin = - 10000;
        float tmax = 10000;
        if (mySphere[k].type == 0) {
            float parameterA = viewDirection.sqrLength();
            float parameterB = 2 * viewDirection.dot(viewOrigin - mySphere[k].center);
            float parameterC = (viewOrigin - mySphere[k].center).sqrLength() - mySphere[k].radius * mySphere[k].radius;
            float parameterD = parameterB * parameterB - 4 * parameterA * parameterC;
            if (parameterD >= 0 && sqrt(parameterD) + parameterB < 0 && ((0 - sqrt(parameterD) - parameterB)/(2 * parameterA) > 0.0001)) {
                t[k] = (0 - sqrt(parameterD) - parameterB)/(2 * parameterA);
            }
            else {
                t[k] = 10000;   //no intersection, set t a large value
            }
        }
        if (mySphere[k].type == 1) {
            float tx1 = (mySphere[k].center[0] - mySphere[k].width * 0.5 - viewOrigin[0])/viewDirection[0];
            float tx2 = (mySphere[k].center[0] + mySphere[k].width * 0.5 - viewOrigin[0])/viewDirection[0];
            float ty1 = (mySphere[k].center[1] - mySphere[k].height * 0.5 - viewOrigin[1])/viewDirection[1];
            float ty2 = (mySphere[k].center[1] + mySphere[k].height * 0.5 - viewOrigin[1])/viewDirection[1];
            float tz1 = (mySphere[k].center[2] - mySphere[k].depth * 0.5 - viewOrigin[2])/viewDirection[2];
            float tz2 = (mySphere[k].center[2] + mySphere[k].depth * 0.5 - viewOrigin[2])/viewDirection[2];
            tmin = max(tmin, min(tx1, tx2));
            tmax = min(tmax, max(tx1, tx2));
            tmin = max(tmin, min(ty1, ty2));
            tmax = min(tmax, max(ty1, ty2));
            tmin = max(tmin, min(tz1, tz2));
            tmax = min(tmax, max(tz1, tz2));
            if ((tmin >= tmax) || (tmin <= 0) ) {
                t[k] = 10000;
            }
            else {
                t[k] = tmin;
            }
            //TODO
        }
    }
    
    float tmin = 9999;
    
    //This is current sphere's number
    int sphere_num = objectNumber + 1; //initial no intersection sphere
    
    for (int k = 0; k < objectNumber; k++) {
        tmin = min(tmin, t[k]);
        if (tmin == t[k]) {
            sphere_num = k; //remember which sphere was intersected
        }
    }
    if (tmin < 9999) {  //generate intersection color
        SbVec3f p = viewOrigin + tmin * viewDirection;
        //mySphere[sphere_num].F.multVecMatrix(p, p);
        SbVec3f Q;
        if (mySphere[sphere_num].type == 0) {
            Q = p - mySphere[sphere_num].center;
            Q.normalize();  // *maybe* This is the intersection's normal vector
        }
        if (mySphere[sphere_num].type == 1) {
            if ((p[0] == mySphere[sphere_num].center[0] - mySphere[sphere_num].width * 0.5) ||(p[0] == mySphere[sphere_num].center[0] + mySphere[sphere_num].width * 0.5)) {
                SbVec3f qx(1, 0, 0);
                Q = qx;
            }
            if ((p[1] == mySphere[sphere_num].center[1] - mySphere[sphere_num].height * 0.5) ||(p[1] == mySphere[sphere_num].center[1] + mySphere[sphere_num].height * 0.5)) {
                SbVec3f qx(0, 1, 0);
                Q = qx;
            }
            if ((p[2] == mySphere[sphere_num].center[2] - mySphere[sphere_num].depth * 0.5) ||(p[2] == mySphere[sphere_num].center[2] + mySphere[sphere_num].depth * 0.5)) {
                SbVec3f qx(0, 0, 1);
                Q = qx;
            }
        }
        
        //try to use solid texturing
        /*SbVec3f blue(0, 0, 1);
        SbVec3f green(0, 1, 0);
        if (sin(PI * p[0] * 4) > 0) {
            return(blue);
        }
        else {
            return(green);
        }*/
        
        bool shadowFlag = true;
        
        //calculate Phong color
        float Ka = 0.2; // assumed by lab2
        SbVec3f phongAmbient(0, 0, 0);
        SbVec3f phongDiffuse(0, 0, 0);
        SbVec3f phongSpecular(0, 0, 0);
        SbVec3f viewVector;    //from intersection to camera;
        SbVec3f *incidentLight;
        SbVec3f *reflectiveLight;
        viewVector = viewOrigin - p;
        viewVector.normalize();
        incidentLight = new SbVec3f[lightNumber];
        reflectiveLight = new SbVec3f[lightNumber];
        
        for (int ii = 0; ii < lightNumber; ii ++) {
            SbVec3f Q1 = Q;
            incidentLight[ii] = myPointLight[ii].location - p;
            incidentLight[ii].normalize();
            reflectiveLight[ii] = 2 * (Q1.dot(incidentLight[ii])) * Q1 - incidentLight[ii];
            phongAmbient += Ka * myPointLight[ii].intensity * mySphere[sphere_num].ambientColor;
            bool intersectFlag = false;
            if (incidentLight[ii].dot(Q1) > 0 && shadowOn == 1) {
                //calculate intersection to decide shadow
                for (int tt = 0; tt < objectNumber; tt ++) {
                    //int tt = 0; (tt < objectNumber && tt != sphere_num); tt++
                    SbVec3f p1 = p; //+ 0.001 * incidentLight[ii];
                    if (mySphere[tt].transparency == 0) { //only non-transparent spheres block light
                        float parameterA = (myPointLight[ii].location - p1).sqrLength();
                        float parameterB = 2 * (myPointLight[ii].location - p1).dot(p1 - mySphere[tt].center);
                        float parameterC = (p1 - mySphere[tt].center).sqrLength() - mySphere[tt].radius * mySphere[tt].radius;
                        float parameterD = parameterB * parameterB - 4 * parameterA * parameterC;
                        if ((parameterD > 0) && ((sqrt(parameterD) + parameterB) < 0)){
                            intersectFlag = true;
                        }
                    }
                }
            }
            if (intersectFlag == false) {
                //(N.L)*diffuseColor*light_intensity*light_color
                //but I didn't use light_color here
                if (sphere_num == (objectNumber - 1)) {
                    if (sin(PI * p[0] * 4) > 0) {
                        phongDiffuse += max(Q1.dot(incidentLight[ii]), (float)0) * mySphere[sphere_num].diffuseColor * myPointLight[ii].intensity;
                        float temp1 = pow(max(viewVector.dot(reflectiveLight[ii]), (float)0), (float)50);
                        phongSpecular += temp1 * mySphere[sphere_num].specularColor * myPointLight[ii].intensity;
                    }
                    else {
                        phongDiffuse += max(Q1.dot(incidentLight[ii]), (float)0) * (white -mySphere[sphere_num].diffuseColor) * myPointLight[ii].intensity;
                        float temp1 = pow(max(viewVector.dot(reflectiveLight[ii]), (float)0), (float)50);
                        phongSpecular += temp1 * (white - mySphere[sphere_num].specularColor) * myPointLight[ii].intensity;
                    }
                }
                else {
                    phongDiffuse += max(Q1.dot(incidentLight[ii]), (float)0) * mySphere[sphere_num].diffuseColor * myPointLight[ii].intensity;
                    float temp1 = pow(max(viewVector.dot(reflectiveLight[ii]), (float)0), (float)50);
                    phongSpecular += temp1 * mySphere[sphere_num].specularColor * myPointLight[ii].intensity;
                }
                //pow((V.R),20)*specularColor*light_intensity*light_color
                //but I didn't use light_color here
                //I use H.N to represent V.R
                //phongSpecular += pow(max(Q1.dot(halfVector[ii]), (float)0), (float)50) * mySphere[sphere_num].specularColor * myPointLight[ii].intensity;
            }
            shadowFlag = !intersectFlag;
        }
        //if (mySphere[sphere_num].transparency == 0) {
        finalColor =  (1 - mySphere[sphere_num].transparency ) * (phongAmbient + phongDiffuse + phongSpecular);
        //}
        //else {
        //    SbVec3f bgcolor(0.3, 0.3, 0.3);
        //    finalColor += bgcolor;
        //}
        //finalColor = (1 - mySphere[sphere_num].transparency) * (phongAmbient + phongDiffuse + phongSpecular);
        //if (mySphere[sphere_num].transparency > 0) {
        //    SbVec3f bgcolor(0.3, 0.3, 0.3);
        //    finalColor += bgcolor;
        //}
        finalColor = scaleColor(finalColor);
        if (recursionDep < MAXRECURSION && shadowFlag) {   // reflectionOn should be here
            SbVec3f D = tmin * viewDirection;
            D.normalize();
            if (mySphere[sphere_num].shininess > 0 && reflectionOn == 1) {   // I move the reflectionOn here
                //SbVec3f D = tmin * viewDirection;
                //D.normalize();
                SbVec3f R = D - 2 * D.dot(Q) * Q;
                R.normalize();
                finalColor += /*mySphere[sphere_num].shininess **/ shade(p, R, recursionDepth + 1);
                finalColor = scaleColor(finalColor);
            }
            
            if (mySphere[sphere_num].transparency > 0) {
                SbVec3f D1 = D;
                D1.negate();    // I want the incident ray point to the camera
                //cout << "!!!!!! " << mySphere[sphere_num].transparency << endl;
                if (D1.dot(Q) > 0) {// that means from air to sphere
                    float index = AIRINDEX / SPHEREINDEX;
                    if (1 - index * index * (1 - (Q.dot(D1))) > 0) {
                        float tempt1 = (index * Q.dot(D1) - sqrt(1 - index * index * (1 - Q.dot(D1) * Q.dot(D1))));
                        SbVec3f T = tempt1 * Q - index * D1;
                        T.normalize();
                        finalColor += mySphere[sphere_num].transparency * shade(p, T, recursionDepth + 1);
                        finalColor = scaleColor(finalColor);
                    }
                    //else {  //that means no refraction 
                    //    //TODO: total internal reflection
                    //}
                }
                else {// that means from sphere to air
                    float index = SPHEREINDEX / AIRINDEX;
                    SbVec3f Qr = mySphere[sphere_num].center - p;
                    Qr.normalize();
                    if (1 - index * index * (1 - (Qr.dot(D))) > 0) {
                        float tempt1 = (index * Qr.dot(D1) - sqrt(1 - index * index * (1 - Qr.dot(D1) * Qr.dot(D1))));
                        SbVec3f T = tempt1 * Qr - index * D1;
                        T.normalize();
                        finalColor += mySphere[sphere_num].transparency * shade(p, T, recursionDepth + 1);
                        finalColor = scaleColor(finalColor);
                    }
                    //TODO
                }
                //finalColor = (1 - mySphere[sphere_num].transparency) * finalColor;
                //finalColor = scaleColor(finalColor);
            }
        }
    }
    else {
        if (recursionDepth == 0) {
            finalColor[0] = 0.3;
            finalColor[1] = 0.3;
            finalColor[2] = 0.3;
        }
        else {
            finalColor[0] = 0;
            finalColor[1] = 0;
            finalColor[2] = 0;
        }
    }
    //finalColor = backgroundColor(finalColor);
    return(finalColor);
}

SbVec3f scaleColor(SbVec3f sourceColor) {
    SbVec3f c = sourceColor;
    float scaleParameter = 1;
    for (int k = 0; k < 3; k ++) {
        scaleParameter = max(c[k], scaleParameter);
    }
    return(c / scaleParameter);
}

SbVec3f backgroundColor(SbVec3f sourceColor) {
    SbVec3f c = sourceColor;
    if(sourceColor[0] == 0 && sourceColor[1] == 0 && sourceColor[2] == 0) {
        sourceColor[0] = 0.3;
        sourceColor[1] = 0.3;
        sourceColor[2] = 0.3;
    }
    return(sourceColor);
}

float randomNo() {
    float i = rand()/(float)(RAND_MAX);
    return(i);
}










