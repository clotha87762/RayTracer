#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include "imageIO.h"
#include "algebra3.h"
#include "hw1.h"

#define epsilon 0.0001

using namespace std;

const float PLANE_DISTANCE = 1.0f;
const float ASPECT_RATIO = 1.0f; // 想要不同的fov就可以設成非0的值  x:y

vec3 upVec(0,1.0f,0);

Camera camera;
vector<SimpleMesh*> simpleMeshes ;



void Sphere::PrintElements(){
        cout<<"center: "<<center[0]<<" "<<center[1]<<" "<<center[2]<<endl;
        cout<<"radius: "<<radius<<endl;
}

bool Sphere::IntersectionTest(vec3 rayDir , vec3 rayOri ){


    vec3 h = center - rayOri;
    float mu = h * rayDir;
    float delta = (mu*mu) - ( h * h) + (radius * radius);

    if(delta < 0){
        //cout<<"sphere false"<<endl;
        return false;
    }
    float t1 = mu + sqrt(delta);
    float t2 = mu - sqrt(delta);

    if(t1<0&&t2<0){
        //cout<<"sphere false"<<endl;
        return false;
    }

    //cout<<"sphere true"<<endl;
    //if( (t1>=0 && t2 < 0) || (t1<0 &&　t2>=0) ) //one intersection, camera in the sphere
     if((t1>=0&&t2<0) || (t1<0 && t2>=0))
        return true;

    return true; //two intersection

}

void Triangle::PrintElements(){
     cout<<"vertices: "<<vertices[0][0]<<" "<<vertices[0][1]<<" "<<vertices[0][2]<<endl;
     cout<<vertices[1][0]<<" "<<vertices[1][1]<<" "<<vertices[1][2]<<endl;
     cout<<vertices[2][0]<<" "<<vertices[2][1]<<" "<<vertices[2][2]<<endl;
}

bool Triangle::IntersectionTest(vec3 rayDir, vec3 rayOri){

    vec3 edge1,edge2 , tvec,pvec,qvec;
    float det,detInv;
    float t,u,v;

    edge1 = vertices[1] - vertices[0];
    edge2 = vertices[2] - vertices[0];

    pvec = rayDir ^ edge2;
    det = edge1 * pvec;


    #ifdef TEST_CULL

        if(det < epsilon){
            return false;
        }
        tvec = rayOri - vertices[0];
        u = tvec * pvec;
        if(u<0||u>det)
            return 0;

        qvec = tvec ^ edge1;

        v = rayDir * qvec;

        if(v<0 || u+v > det)
            return 0;

        t = edge2 ^ qvec;
        detInv = 1.0f / det;
        t *= detInv;
        u *= detInv;
        v *= detInv;

    #else

        if(det > -epsilon && det < epsilon)
            return false;
        detInv = 1.0f / det;
        tvec = rayOri - vertices[0];

        u = (tvec*pvec) * detInv;
        if(u<0||u>1.0f)
            return 0;

        qvec = tvec ^ edge1;
        v = (rayDir * qvec) * detInv;
        if(v < 0 || u+v > 1.0f)
            return false;

        t = (edge2 * qvec) * detInv;

    #endif // TEST_CULL

        return true;
}

void Camera::PrintElements(){
            cout<<center[0]<<" "<<center[1]<<" "<<center[2]<<endl;
            cout<<viewDirection[0]<<" "<<viewDirection[1]<<" "<<viewDirection[2]<<endl;
            cout<<fovX<<"  "<<fovY<<endl;
            cout<<planeDistance<<endl;
            cout<<resX<<endl;
            cout<<resY<<endl;
}

void Camera::CalcPlaneParameters(){

    planeCenter = center + ((viewDirection/viewDirection.length()) * planeDistance);
    planeXAxis  = (up ^ viewDirection).normalize();
    planeYAxis = (planeXAxis ^ viewDirection).normalize();
    planeNormal = viewDirection;

    cout<<up[0]<<" "<<up[1]<<" "<<up[2]<<endl;
    cout<<viewDirection[0]<<" "<<viewDirection[1]<<" "<<viewDirection[2]<<endl;
    cout<<planeXAxis[0]<<" "<<planeXAxis[1]<<" "<<planeXAxis[2]<<endl;
    cout<<planeYAxis[0]<<" "<<planeYAxis[1]<<" "<<planeYAxis[2]<<endl;

    float thetaX = fovX / 2.0f;
    float thetaY = fovY / 2.0f;
    float halfWidth = planeDistance * tan( thetaX * (3.14159f/180.0f));
    float halfHeight = planeDistance * tan( thetaY * (3.14159/180.0f) );

    cout<<halfWidth<<"  "<<halfHeight<<endl;

    planeCorners[0] =  planeCenter -  (planeXAxis * halfWidth) - (planeYAxis*halfHeight);
    planeCorners[1] =  planeCenter +  (planeXAxis * halfWidth) - (planeYAxis*halfHeight);
    planeCorners[2] =  planeCenter -  (planeXAxis * halfWidth) + (planeYAxis*halfHeight);
    planeCorners[3] =  planeCenter +  (planeXAxis * halfWidth) + (planeYAxis*halfHeight);

}

vec3 Camera::GetPixelCenter(vec2 offset){
    float pixelWidth = (planeCorners[1] - planeCorners[0]).length() / ((float)resX);
    float pixelHeight =  (planeCorners[0] - planeCorners[2]).length() / ((float)resY);

    //cout<<"pixelWidth "<<pixelWidth<<" height "<<pixelHeight<<endl;

    return planeCorners[0] + (0.5f * pixelWidth * planeXAxis) + (0.5f * pixelHeight * planeYAxis) +
    (offset[0] * pixelWidth * planeXAxis) + (offset[1] * pixelHeight * planeYAxis);
}

int main(int argc , char* argv[]){

    char* inputName = "input.txt";
    char* str = new char[100];


    float f;
    int i;


    if(argc>1){
        inputName = argv[1];
    }


    ifstream input (inputName,ifstream::in);

    if(!input.is_open()){
        cout<<"can not open input file"<<endl;
        return -1;
    }


    //camera = new Camera();
    camera.planeDistance = PLANE_DISTANCE;
    camera.up = upVec;

    try{


        while(input >> str){
            cout<<str<<endl;
            if( !strcmp(str,"E")){
                //camera.center = new vec3();
                input >> str;
                camera.center[0] = atof(str);
                input >> str;
                camera.center[1] = atof(str);
                input >> str;
                camera.center[2] = atof(str);

            }
            else if(!strcmp(str ,"V")){
                //camera.viewDirection = new vec3();

                input >> str;
                camera.viewDirection[0] = atof(str);
                input >> str;
                camera.viewDirection[1]  = atof(str);
                input >> str;
                camera.viewDirection[2]  = atof(str);
                camera.viewDirection= camera.viewDirection.normalize();

            }
            else if( !strcmp(str ,"F")){
                input >> str;
                f = (float)atof(str);
                camera.fovX = f;
                camera.fovY = f/(ASPECT_RATIO);
            }
            else if(!strcmp(str,"R")){
                input >> str;
                i = atoi(str);
                camera.resX = i;
                input >> str;
                i = atoi(str);
                camera.resY = i;
            }
            else if(!strcmp(str ,"S")){
                Sphere* sphere = new Sphere();
                //sphere.center = new vec3();
                input >> str;
                sphere->center[0] = atof(str);
                input >> str;
                sphere->center[1] = atof(str);
                input >> str;
                sphere->center[2] = atof(str);
                input >> str;
                sphere->radius = atof(str);
                simpleMeshes.push_back(sphere);
            }
            else if( !strcmp(str,"T")){
                Triangle* triangle = new Triangle();
                //triangle.vertices[0] = new vec3();
                //triangle.vertices[1] = new vec3();
                //triangle.vertices[2] = new vec3();
                input >> str;
                triangle->vertices[0][0] = atof(str);
                input >> str;
                triangle->vertices[0][1] = atof(str);
                input >> str;
                triangle->vertices[0][2] = atof(str);
                input >> str;
                triangle->vertices[1][0] = atof(str);
                input >> str;
                triangle->vertices[1][1] = atof(str);
                input >> str;
                triangle->vertices[1][2] = atof(str);
                input >> str;
                triangle->vertices[2][0] = atof(str);
                input >> str;
                triangle->vertices[2][1] = atof(str);
                input >> str;
                triangle->vertices[2][2] = atof(str);
                simpleMeshes.push_back(triangle);

            }

        }
    }
    catch(...){
        cout<< "Wrong input file format!" <<endl;
    }

    camera.CalcPlaneParameters();
    GenerateBinaryIntersectionTestImage(camera);

    /*
    camera.PrintElements();
    for(i=0;i<simpleMeshes.size();i++){
        simpleMeshes[i]->PrintElements();
    }
    */

     for(i=0;i<simpleMeshes.size();i++){
        delete simpleMeshes[i];
    }

    delete[] str;
	return 0;
}


bool RayCastInterSectionTestSimpleMesh(Camera  cam , vec2 xy ){
    vec3 pixelCenter = cam.GetPixelCenter(xy);
    vec3 ray = (pixelCenter-cam.center).normalize();

    //cout<<"p "<<pixelCenter[0]<<" "<<pixelCenter[1]<<" "<<pixelCenter[2]<<endl;
    //cout<<"c "<<cam.center[0]<<" "<<cam.center[1]<<" "<<cam.center[2]<<endl;

    vector<SimpleMesh*>::iterator meshIte;
    bool b = false;

    for(meshIte=simpleMeshes.begin(); meshIte!=simpleMeshes.end(); ++meshIte){

        b = (*meshIte)->IntersectionTest(ray,cam.center);

        if(b) break;
    }

    return b;
}


void GenerateBinaryIntersectionTestImage(Camera cam){

	ColorImage image;
	int x, y;
	Pixel p={0,0,0};

	image.init(cam.resX, cam.resY);
    vec2 vv;
    int q1 = 0 ,q2=0;
	for(x =0;x<cam.resX;x++){
        for(y=0;y<cam.resY;y++){
            vv[0] = x;
            vv[1] = y;
            if(RayCastInterSectionTestSimpleMesh(cam,vv)){
                //cout<<"draw white"<<q1++<<endl;
                p.R = 255;
                p.G = 255;
                p.B = 255;
            }
            else{
                //cout<<"draw black"<<q2++<<endl;
                p.R =0;
                p.G =0;
                p.B =0;
            }
            image.writePixel(x,y,p);
        }
	}
    /*
	for (y=0; y<cam.resY; y++) {
		for (x=0; x<256; x++) {
			p.R = y;
			image.writePixel(x, y, p);
		}
	}
	*/

	char s[20] = "output.ppm";
	image.outputPPM(s);

}
