#ifndef __Camera_h__
#define __Camera_h__

#include "algebra3.h"

class Camera{

public :

    int ssaaRatio;

    vec3 center;
    vec3 up;
    vec3 viewDirection;
    float fovX , fovY ;
    float planeDistance;
    int resX;
    int resY;

    vec3 planeCenter;
    vec3 planeNormal;
    vec3 planeXAxis;
    vec3 planeYAxis;
    vec3 planeCorners[4]; // topLeft ,topRight, bottomLeft , bottomRight


    Camera(){
        ssaaRatio = 1;
        planeDistance = 1.0f; // planeDistance好像沒在input格式中定義...QAQ
    }

    Camera(vec3 center , vec3 up, vec3 viewDirection,float fovX,float fovY,int resX,int resY){
        this->center = center;
        this->up = up;
        this->viewDirection = viewDirection;
        this->fovX = fovX;
        this->fovY = fovY;
        this->resX = resX;
        this->resY = resY;
    }

    vec3 GetPixelCenter(vec2 offset);
    void PrintElements();
    void CalcPlaneParameters();
};
#endif // __Camera_h__
