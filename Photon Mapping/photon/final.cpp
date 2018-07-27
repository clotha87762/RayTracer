#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <ctime>
#include <stack>
#include <algorithm>


#include "imageIO.h"
#include "algebra3.h"
#include "photon.h"
#include "parameter.h"
#include "final.h"




using namespace std;

// =========== Adjustable parameters ===================
double totalIntersectionCheck = 0;
double averageCheck = 0;
int checkTime = 0;

double elapsedTime = 0;

//KDNode* rootKD;
Camera camera;
vector<SimpleMesh*> simpleMeshes ;
vector<SimpleLight*> simpleLights;

vector<SimpleMesh*> causticMeshes;

vector<SimpleObject*> simpleObjects;

vector<SimpleMesh*> lightMeshes;

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

/*
    if(t1<0||t2<0){
        //cout<<"sphere false"<<endl;
        return false;
    }
*/
    if(t1<=0 && t2 <=0)
        return false;
    //cout<<"sphere true"<<endl;

    if( t1>=0 && t2 < 0){
        tempIntersectionPoint = rayOri + t1*rayDir;
        return true;
    }else if(t1<0&&t2>=0){
        tempIntersectionPoint = rayOri + t2*rayDir;
        return true;
    } //one intersection, camera in the sphere
    else {

        float t = t1 < t2? t1:t2;
        tempIntersectionPoint = rayOri + t*rayDir;

    }

    float t = t1 < t2? t1:t2;
    tempIntersectionPoint = rayOri + t*rayDir;

    return true; //two intersection

}

vec3 Sphere::GetNormal(vec3 point){
    return (point - center).normalize();
}

vec3 Triangle::GetNormal(vec3 point){

    if(BARYCENTRIC_INTERPOLATION_NORMAL){
        float a0,a1,a2,A;
        a0 = ((vertices[1] - point) ^ (vertices[2] - point)).length();
        a1 = ((vertices[0] - point) ^ (vertices[2] - point)).length();
        a2 = ((vertices[0] - point) ^ (vertices[1] - point)).length();
        A = ((vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0])).length();
        vec3 v = ((a0/A) * vertices[0]) + ((a1/A) * vertices[1])  + ((a2/A) * vertices[2]) ;
        return v.normalize();
    }
    else{

        return normal[0].normalize(); // 不做內差則直接回傳
    }
}

vec3 Sphere::GetIntersectionPoint(vec3 rayDir, vec3 rayOri){


}

vector<vec3> Sphere::GetPositions(){
    vector<vec3> toReturn;
    toReturn.push_back(center);

    return toReturn;
}

vec3 Triangle::GetIntersectionPoint(vec3 rayDir, vec3 rayOri){


}


vector<vec3> Triangle::GetPositions(){
    vector<vec3> toReturn;
    toReturn.push_back(vertices[0]);
    toReturn.push_back(vertices[1]);
    toReturn.push_back(vertices[2]);
    return toReturn;
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
            return false;

        qvec = tvec ^ edge1;
        v = (rayDir * qvec) * detInv;
        if(v < 0 || u+v > 1.0f)
            return false;

        t = (edge2 * qvec) * detInv;
        if(t<0)
            return false;

    #endif // TEST_CULL

        tempIntersectionPoint = rayOri + (rayDir*t);
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


int KDNode::index = 0;

bool KDNode::HitTest(vec3 ori , vec3 ray)
{
   // if(kdMeshes.size()<=0)
     //   return false;
    //cout<<"test"<<endl;

    float tmin = ((center[0]-extend[0]/2.0f) - ori[0]) / ray[0];
    float tmax = ((center[0]+extend[0]/2.0f) - ori[0]) / ray[0];

    if (tmin > tmax) swap(tmin, tmax);

    float tymin = ((center[1]-extend[1]/2.0f) - ori[1]) / ray[1];
    float tymax = ((center[1]+extend[1]/2.0f) - ori[1]) / ray[1];

    if (tymin > tymax) swap(tymin, tymax);

    if ((tmin > tymax) || (tymin > tmax)){
       //cout<<"not hit: ray "<<ray[0]<<" "<<ray[1]<<" "<<ray[2]<<" ori"<<ori[0]<<" "<<ori[1]<<" "<<ori[2]<<endl;
       // cout<<"center "<<center[0]<<" "<<center[1]<<" "<<center[2]<<" ext "<<extend[0]<<" "<<extend[1]<<" "<<extend[2]<<endl;
        return false;
    }

    if (tymin > tmin)
    tmin = tymin;

    if (tymax < tmax)
    tmax = tymax;

    float tzmin = ((center[2]-extend[2]/2.0f) - ori[2]) / ray[2];
    float tzmax = ((center[2]+extend[2]/2.0f) - ori[2]) / ray[2];

    if (tzmin > tzmax) swap(tzmin, tzmax);

    if ((tmin > tzmax) || (tzmin > tmax)){
        //cout<<"not hit: ray "<<ray[0]<<" "<<ray[1]<<" "<<ray[2]<<" ori"<<ori[0]<<" "<<ori[1]<<" "<<ori[2]<<endl;
        //cout<<"center "<<center[0]<<" "<<center[1]<<" "<<center[2]<<" ext "<<extend[0]<<" "<<extend[1]<<" "<<extend[2]<<endl;
        return false;
    }

    if (tzmin > tmin)
    tmin = tzmin;

    if (tzmax < tmax)
    tmax = tzmax;

    //cout<<"hitt"<<endl;
    return true;
}

void KDNode::CalcBound(){
    vec3 mins(1e9,1e9,1e9);
    vec3 maxs(-1e9,-1e9,-1e9);


    for(int i=0;i<kdMeshes.size();i++){

        Triangle* tri = dynamic_cast<Triangle*>(kdMeshes[i]);
        if(tri != NULL){
            for(int j=0;j<3;j++){
               for(int k=0;k<3;k++){
                    if(tri->vertices[j][k] < mins[k]){
                        mins[k] = tri->vertices[j][k];
                    }
                    if(tri->vertices[j][k] > maxs[k]){
                        maxs[k] = tri->vertices[j][k];
                    }

                }
            }
        }

        Sphere* sph = dynamic_cast<Sphere*>(kdMeshes[i]);
        if(sph != NULL){
            for(int j=0;j<3;j++){
                if(sph->center[j] + sph->radius > maxs[j]){
                    maxs[j] = sph->center[j] + sph->radius;
                }
                if(sph->center[j] - sph->radius < mins[j]){
                    mins[j] = sph->center[j] - sph->radius;
                }
            }
        }



    }

     center = (mins + maxs)/2.0f;
     extend = maxs-mins;

     /*
     for(int i=0;i<3;i++)
     if(extend[i]<boundingBoxMinimunSize){
            extend[i] = boundingBoxMinimunSize;

     }
     */

    //cout<<"center "<<center[0]<<" "<<center[1]<<" "<<center[2]<<" ext"<<extend[0]<<" "<<extend[1]<<" "<<extend[2]<<endl;

    float temp = -1e9;
    for(int i=0;i<3;i++){
        if( extend[i] > temp){
            temp = extend[i];
            splitAxis = i;
        }
    }
    //cout<<"split axis "<<splitAxis<<endl;
}

int builded = 0;
int maxDepth = 0;
void KDNode::BuildKD(int depth){

    layer = depth;
    if(depth > maxDepth)
        maxDepth = depth;

    //if(depth>=kdMaxDepth)
    //    return;
    cout<<"builded "<<builded++<<endl;

    CalcBound();

    vector<float> poses;
    for(int i=0;i<kdMeshes.size();i++){
        Triangle* tri = dynamic_cast<Triangle*>(kdMeshes[i]);
        Sphere* sph = dynamic_cast<Sphere*>(kdMeshes[i]);
        if(tri !=NULL){
            poses.push_back(tri->vertices[0][splitAxis]);
            poses.push_back(tri->vertices[1][splitAxis]);
            poses.push_back(tri->vertices[2][splitAxis]);
        }
        if(sph !=NULL){
            poses.push_back(sph->center[splitAxis]);
        }
    }
    std::sort(poses.begin(),poses.end());
    splitValue = poses[poses.size()/2];

    //cout<<"split value:"<<splitValue<<endl;

    leftChild = new KDNode();
    rightChild = new KDNode();
    //if(depth<kdMaxDepth){
        int match = 0;
        for(int i=0;i<kdMeshes.size();i++){

            Triangle* tri = dynamic_cast<Triangle*>(kdMeshes[i]);
            Sphere* sph = dynamic_cast<Sphere*>(kdMeshes[i]);
            if(tri !=NULL){
               bool inleft,inright;
                inleft = false;
                inright = false;
               for(int j=0;j<3;j++){
                if(tri->vertices[j][splitAxis] >= splitValue){
                    rightChild->kdMeshes.push_back(tri);
                    rightChild->nonZero = true;
                    inright = true;
                    break;
                }
               }
               for(int j=0;j<3;j++){
                if(tri->vertices[j][splitAxis] < splitValue){
                    leftChild->kdMeshes.push_back(tri);
                    leftChild->nonZero = true;
                    inleft = true;
                    break;
                }
               }
                if(inleft&&inright)
                    match++;
            }
            if(sph !=NULL){
                cout<<"sphere??"<<endl;
                if(sph->center[splitAxis] + sph->radius >= splitValue){
                    rightChild->kdMeshes.push_back(sph);
                }
                if(sph->center[splitAxis] - sph->radius < splitValue){
                    leftChild->kdMeshes.push_back(tri);
                }
            }
        }

        //cout<<"match ratio:"<<(float)match/(float)leftChild->kdMeshes.size()<<endl;
        //cout<<"match "<<match<<" left size"<<leftChild->kdMeshes.size() <<" right size"<<rightChild->kdMeshes.size()<<endl;


        if((float)match<(float)leftChild->kdMeshes.size()/modelScale && (float)match<(float)rightChild->kdMeshes.size()/modelScale ){

                leftChild->BuildKD(depth+1);
                rightChild->BuildKD(depth+1);
                splitLeft = true;
                splitRight = true;

        }
        else{
            leftChild->CalcBound();
            rightChild->CalcBound();
            leftChild->layer = depth+1;
            rightChild->layer = depth+1;
            //leftChild->isLeaf = true;
            //rightChild->isLeaf = true;
            this->isLeaf = true;
            //isLeaf = true;
            cout<<"max depth"<<depth+1<<endl;
        }
    //}

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

    bool materialIncluded = false;
    vec3 tempSurfColor(1.0f,1.0f,1.0f);
    float tempKs =1.0f, tempKd =1.0f , tempKa =1.0f;
    float tempSpecularity = 6.0f;
    float tempReflect = 0.9f;
    float tempRefract = 0;
    float tempNr = 1.0f;
    //float tempExp = 6.0f;

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
                camera.resX = i *SSAA_RATIO ;
                input >> str;
                i = atoi(str);
                camera.resY = i * SSAA_RATIO;

                camera.ssaaRatio = SSAA_RATIO;

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

                sphere->Ka = tempKa;
                sphere->Kd = tempKd;
                sphere->Ks = tempKs;
                sphere->surfColor = tempSurfColor;
                sphere->specularity = tempSpecularity;
                sphere->reflect = tempReflect;
                sphere->refract = tempRefract;
                sphere->Nr = tempNr;

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

                triangle->Ka = tempKa;
                triangle->Kd = tempKd;
                triangle->Ks = tempKs;
                triangle->surfColor = tempSurfColor;
                triangle->specularity = tempSpecularity;
                triangle->reflect = tempReflect;
                triangle->refract = tempRefract;
                triangle->Nr = tempNr;



                vec3 tempCross = (triangle->vertices[2] -triangle->vertices[0])^ (triangle->vertices[1] -triangle->vertices[0]);
                triangle->normal[0] = tempCross.normalize();
                triangle->normal[1] = tempCross.normalize();
                triangle->normal[2] = tempCross.normalize();
                cout<<"n:"<<triangle->normal[0][0]<<" "<<triangle->normal[0][1]<<" "<<triangle->normal[0][2]<<endl;

                simpleMeshes.push_back(triangle);

            }
            else if(!strcmp(str,"L")){
                 PointLight* pl = new PointLight();
                 input >> str;
                 pl->pos[0] = atof(str);
                 input >> str;
                 pl->pos[1] = atof(str);
                 input >> str;
                 pl->pos[2] = atof(str);

                 simpleLights.push_back(pl);
            }
            else if(!strcmp(str,"CL")){  //Color Light

            }
            else if(!strcmp(str,"DL")){

            }
            else if(!strcmp(str,"SL")){

            }
            else if(!strcmp(str,"ML")){  // ML : 有加Color和intensity資訊
                 PointLight* pl = new PointLight();
                 input >> str;
                 pl->pos[0] = atof(str);
                 input >> str;
                 pl->pos[1] = atof(str);
                 input >> str;
                 pl->pos[2] = atof(str);
                 input >> str;
                 pl->color[0] = atof(str);
                 input >> str;
                 pl->color[1] = atof(str);
                 input >> str;
                 pl->color[2] = atof(str);
                 input >> str;
                 pl->intensity = atof(str);
                 simpleLights.push_back(pl);
            }
            else if(!strcmp(str,"MDL")){

            }
            else if(!strcmp(str,"MSL")){

            }
            else if(!strcmp(str,"M")){
                 input >> str;
                 tempSurfColor[0] = atof(str);
                 input >> str;
                 tempSurfColor[1] = atof(str);
                 input >> str;
                 tempSurfColor[2] = atof(str);
                 input >> str;
                 tempKa = atof(str);
                 input >> str;
                 tempKd = atof(str);
                 input >> str;
                 tempKs = atof(str);
                 input >> str;
                 tempSpecularity = atof(str);
                 input >> str;
                 tempReflect = atof(str);
                 input >> str;
                 tempRefract = atof(str);
                input >> str;
                 tempNr = atof(str);
            }
            else if(!strcmp(str,"RectL")){

                 RectLight* rl = new RectLight();
                 input >> str;
                 rl->pos[0] = atof(str);
                 input >> str;
                 rl->pos[1] = atof(str);
                 input >> str;
                 rl->pos[2] = atof(str);

                 rl->dir1 = vec3( 1.0,0,0);
                 rl->dir2 = vec3( 0 , 0 , 1.0);

                 input >> str;
                 rl->size1 = atof(str);
                 input >> str;
                 rl->size2 = atof(str);

                 input >> str;
                 rl->sampleCountDir1 = atoi(str);

                 input >> str;
                 rl->sampleCountDir2 = atoi(str);

                 simpleLights.push_back(rl);

                Triangle* triangle = new Triangle();

                triangle->vertices[0] = rl->pos - rl->dir1*rl->size1*0.5 - rl->dir2*rl->size2*0.5;
                triangle->vertices[1] = rl->pos + rl->dir1*rl->size1*0.5 - rl->dir2*rl->size2*0.5;
                triangle->vertices[2] = rl->pos - rl->dir1*rl->size1*0.5 + rl->dir2*rl->size2*0.5;

                vec3 tempCross = (triangle->vertices[2] -triangle->vertices[0])^ (triangle->vertices[1] -triangle->vertices[0]);
                triangle->normal[0] = tempCross.normalize();
                triangle->normal[1] = tempCross.normalize();
                triangle->normal[2] = tempCross.normalize();
                lightMeshes.push_back(triangle);

                triangle = new Triangle();

                triangle->vertices[0] = rl->pos + rl->dir1*rl->size1*0.5 + rl->dir2*rl->size2*0.5;
                triangle->vertices[1] = rl->pos + rl->dir1*rl->size1*0.5 - rl->dir2*rl->size2*0.5;
                triangle->vertices[2] = rl->pos - rl->dir1*rl->size1*0.5 + rl->dir2*rl->size2*0.5;

                tempCross = (triangle->vertices[2] -triangle->vertices[0])^ (triangle->vertices[1] -triangle->vertices[0]);
                triangle->normal[0] = tempCross.normalize();
                triangle->normal[1] = tempCross.normalize();
                triangle->normal[2] = tempCross.normalize();
                lightMeshes.push_back(triangle);

            }

        }
    }
    catch(...){
        cout<< "Wrong input file format!" <<endl;
    }

    camera.CalcPlaneParameters();

/*
    for(int i=0;i<simpleMeshes.size();i++){
        simpleMeshes[i]->PreComputeAmbientDiffuse();
    }
*/
    //GenerateBinaryIntersectionTestImage(camera);
    //GenerateColorImage(camera);


    /* 設定ｃａｕｓｔｉｃ區域）*/

    Sphere* causticSphere = new Sphere();
    causticSphere->center = vec3(0.2, -0.3, 0.45);
    causticSphere->radius = 0.2;
    //causticSphere->
    causticMeshes.push_back(causticSphere);


    InitDirectionTable();
    PhotonMappingImageSynthesis(camera);
    /*
    camera.PrintElements();
    for(i=0;i<simpleMeshes.size();i++){
        simpleMeshes[i]->PrintElements();
    }
    */

     for(i=0;i<simpleMeshes.size();i++){
        delete simpleMeshes[i];
    }

    for(i=0;i<causticMeshes.size();i++){
        delete causticMeshes[i];
    }

    delete[] str;
	return 0;
}




vec3 SimpleMeshRayTrace(PhotonMap& pmap,PhotonMap& cpmap , vec3 ori , vec3 dir , int depth ,Camera cam ,float currentNr ){


    //cout<<"depth "<<depth<<endl;
    if(depth>5)
        return vec3(0,0,0);

    float nearestDis = 1e9;

    SimpleMesh* nearestMesh = 0;
    vector<SimpleMesh*>::iterator meshIte;

    vec3 colour(0,0,0);
    vec3 hitPoint;

    checkTime++;

    for(int i=0;i<simpleMeshes.size();i++){

         bool b = simpleMeshes[i]->IntersectionTest(dir,ori);
         if(b){
            vec3 v = simpleMeshes[i]->tempIntersectionPoint;
            float len = (ori-v).length();
            if(len<nearestDis && len > TOLERANCE_DEPTH && (culling==false || simpleMeshes[i]->GetNormal(v)*dir<0)){
                //cout<<"hit  dis:"<<len<<endl;
                nearestMesh = simpleMeshes[i];
                nearestDis = len;
                hitPoint = v;

            }
         }
    }


    //cout<<"======"<<endl;

    vec3 ambient(0,0,0);
    vec3 diffuse(0,0,0);
    vec3 specular(0,0,0);



    if(nearestMesh){

        if(depth==0)
            rayHitAnything = true;

        //nearestMesh->PrintElements();

        vec3 normal = nearestMesh->GetNormal(hitPoint);
        vec3 normalL = normal;
        if(normal*dir > 0 ){
            normalL = normal* -1.0f;
        }

        vec3 reflected;
        vec3 transmitted;

          if(nearestMesh->reflect > 0.0001f){

            float dot = (dir*normalL);
            reflected = dir - 2 * dot * normalL;

            vec3 reflectColor = /*nearestMesh->surfColor **/ SimpleMeshRayTrace(pmap,cpmap,hitPoint,reflected ,depth+1,cam ,currentNr ); // reflected是否要反過來???
            colour = colour + reflectColor * nearestMesh->reflect;
            //cout<<"ref col "<<colour[0]<<" "<<colour[1]<<" "<<colour[2]<<endl;
        }


        if(nearestMesh->refract > 0.0001){

            /*
            float n1 = currentNr;
            float n2 = nearestMesh->Nr;
            float cos1 = dir * normal;
            float sin2 = (n1/n2) * sqrt( 1-(cos1*cos1));
            transmitted = (   ((n1/n2)*dir) + ( cos1*(n1/n2)  - sqrt(1-(sin2*sin2)))*normal   ).normalize();
            //cout<<transmitted*dir<<endl;

            vec3 displacePoint = hitPoint - normal * refractBias;

            vec3 refractColor = SimpleMeshRayTrace(pmap, cpmap, displacePoint, transmitted ,depth+1,cam , n2 ); // reflected是否要反過來???
             //cout<<refractColor[0]<<"  "<<refractColor[1]<<"  "<<refractColor[2]<<endl;
             colour = colour + refractColor * nearestMesh->refract;
            //cout<<"refracttt"<<endl;
            */
            vec3 reflRay(dir-normal*2*normal*dir);     // Ideal dielectric REFRACTION
            bool into = normal * normalL >0;                // Ray from outside going in?
            float nc=1, nt=currentNr, nnt=into?nc/nt:nt/nc, ddn=dir * normalL, cos2t;

            if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    {// Total internal reflection
               //cout<<"what"<<endl;
                vec3 refractColor =   SimpleMeshRayTrace(pmap,cpmap,hitPoint,reflRay ,depth+1,cam ,currentNr );
                colour = colour + refractColor * nearestMesh->refract;
                //f.mult(raytrace(pmap, pmap_caustic, origin_new,reflRay,triangle_list,light_list,depth,RI));
            }
            else{
                //cout<<"transmit"<<endl;
                vec3 tdir = (dir*nnt - normal*((into?1:-1)*(ddn*nnt+sqrt(cos2t))));
                tdir = tdir.normalize();
                float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir*normal);
                float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re;


                //Re = 0;
                //Tr = 0;
                //cout<<"surf color "<<nearestMesh->surfColor[0]<<" "<<nearestMesh->surfColor[1]<<" "<<nearestMesh->surfColor[2]<<endl;
                //cout<<"Re "<<Re<<"  Tr "<<Tr<<endl;
                vec3 refractColor = /*nearestMesh->surfColor * */
                (SimpleMeshRayTrace(pmap,cpmap,hitPoint,reflRay ,depth+1,cam ,currentNr )*Re+
                 SimpleMeshRayTrace(pmap,cpmap,hitPoint,tdir ,depth+1,cam ,currentNr )*Tr
                );
                colour = colour + refractColor * nearestMesh->refract ;
            }
        }

        //ambient = (nearestMesh->surfColor * ambientIntensity) * nearestMesh->Ka;
        //colour = colour+ambient;

        /* Caustic effect */

        //float color[3];
		//float pos[3]={origin_new.x,origin_new.y,origin_new.z};
		//float normal[3] = {n.x, n.y, n.z};
		vec3 caustic = cpmap.CalcIrradiance( hitPoint, normal , 0.1,  PHOTON_QUERY_COUNT , MIN_PHOTON_TO_ESTIMATE );
        colour = colour + caustic ;

        /* Global illumination */



        for(int i=0;i<GLOBAL_ILLUMINATION_SAMPLE;i++){
            //cout<<"sample "<<i<<endl;
            double r1=2.0*PI*mRandom(), r2=mRandom(), r2s=sqrt(r2);

     		vec3  w=normalL;
            vec3  u=( (fabs(w[0])>0.1 ? vec3(0,1.0f,0):vec3(1.0f,0,0))^w);
            u = u.normalize();
            vec3 v=w^u;

            vec3 d = (u*cos(r1)*r2s + v*sin(r1)*r2s +  w*sqrt(1-r2));
            d = d.normalize();

            vec3 temp = TraceIrradiance(pmap,simpleMeshes , hitPoint, d, depth+1 , currentNr)/(float)GLOBAL_ILLUMINATION_SAMPLE;
            //cout<<"temp "<<temp[0]<<"  "<<temp[1]<<"  "<<temp[2]<<endl;
            colour = colour + temp;
        }




        /*  Direct illumination */

        //ambient = (nearestMesh->surfColor * ambientIntensity) * nearestMesh->Ka;
        //colour = colour+ambient;




        for(int i=0;i<simpleLights.size();i++){

            vector<vec3> samplePoses = simpleLights[i]->GetSamplePoses();
            float factor = 1.0 / (samplePoses.size()>=1 ?samplePoses.size():1);
            //cout <<"factor "<<factor<<"  sample size "<<samplePoses.size()<<endl;

            for(int j=0;j<samplePoses.size();j++){

                //cout<<"sample pos "<<samplePoses[j][0]<<" "<<samplePoses[j][1]<<" "<<samplePoses[j][2]<<endl;

                if(  !IsInShadow(hitPoint , samplePoses[j] , nearestMesh ,normalL) ){
                       //cout<<"qq2"<<endl;
                       vec3 toLight = (samplePoses[j] - hitPoint).normalize();

                       float dot = max(toLight*normalL,0.0f);
                       diffuse = nearestMesh-> Kd * dot *(simpleLights[i]->intensity*nearestMesh->surfColor) ;

                       //cout<<"dot "<<dot<<endl;

                       vec3 eyeVec = (cam.center - hitPoint).normalize();
                       vec3 halfVec = (eyeVec + toLight).normalize();
                       dot = max(halfVec*normalL,0.0f);
                       specular = nearestMesh->Ks * (simpleLights[i]->intensity*nearestMesh->surfColor) * pow((dot),nearestMesh->specularity);

                }else{

                }
                 //colour = colour + ambient + diffuse + specular;
                 colour = colour + factor * (diffuse + specular);
            }
            //cout<<"direct color "<<colour[0]<<" "<<colour[1]<<" "<<colour[2]<<endl;
        }



    }


    return colour;
}

bool IsInShadow(vec3 point , vec3 lightPos , SimpleMesh* currentMesh , vec3 pointNormal){



    point = point + pointNormal * shadowBias ;

    vec3 rayDir = (lightPos - point).normalize();
    float nearestDis = 1e9;
    vector<SimpleMesh*>::iterator meshIte;

    float tLight = (lightPos - point).length();

    for(meshIte=simpleMeshes.begin(); meshIte!=simpleMeshes.end(); ++meshIte){

         if((*meshIte)->refract>0.1)
            continue;

         bool b;
         b = (*meshIte)->IntersectionTest(rayDir,point);
         if(b){
            vec3 v = (*meshIte)->tempIntersectionPoint;
            float len = (point-v).length();
            if(len< tLight && len> TOLERANCE_DEPTH ){ //避免下一次ray trace 馬上又跟自己的mesh intersect到
                return true;
            }
         }

    }

    checkTime++;

    return false;

}


bool TestVisualizeLight(vec3 ori , vec3 dir){


    for(int i=0;i<lightMeshes.size();i++){
        if(lightMeshes[i]->IntersectionTest(dir,ori)){
            //cout<<"yoooooo"<<endl;
            return true;
        }
    }
    return false;
}

float AAFilter[KERNEL_SIZE][3] 	=		/* X, Y, coef */
{

		-0.52, 0.38, 0.128,
		0.41, 0.56, 0.119,
		0.27, 0.08, 0.294,
		-0.17, -0.29, 0.249,
		0.58, -0.55, 0.104,
		-0.31, -0.71, 0.106

};

void PhotonMappingImageSynthesis(Camera cam){

    int numPhoton = 100000;
    int numCausticPhoton = 1000;

    PhotonMap pmap(numPhoton);
    PhotonMap cpmap(numCausticPhoton);

    for(int i=0;i<simpleLights.size();i++){

        //燈一定要放在天花板上

        while (pmap.storedPhotonNum < (float) pmap.maxPhotons / (float) simpleLights.size()  ) {


			float offset_x = mRandom();
			float offset_y = mRandom();

			vec3 pos = simpleLights[i]->pos  - simpleLights[i]->size1*simpleLights[i]->dir1 *0.5+  simpleLights[i]->size1*simpleLights[i]->dir1*offset_x
								- simpleLights[i]->size2*simpleLights[i]->dir2*0.5   +  simpleLights[i]->size2 * simpleLights[i]->dir2*offset_y;
			//vec3 pos = simpleLights[i]->pos;

			vec3 nl(0, -1.0 ,0);

			float r1= 2*PI*mRandom(), r2=mRandom(), r2s=sqrt(r2);

			vec3 w=nl, u= vec3(1.0,0,0)^w ;
			u = u.normalize();
			vec3 v =w^u;

			vec3 d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2));
			d = d.normalize();

			vec3 color(3.0,3.0,3.0);

			PhotonTrace( pmap, simpleMeshes , pos , d , 0 , 2.40 , color );
            //cout<<"stored photon "<<pmap.storedPhotonNum<<endl;
            //cout<<"pos "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
            //cout<<"d "<<d[0]<<" "<<d[1]<<" "<<d[2]<<endl;

			//photon_tracing(pmap_r,pos,d,surfaces,color,0,2.40);
		}

         cout<<"finish generating photons"<<endl;

         while (cpmap.storedPhotonNum < (float) cpmap.maxPhotons / (float) simpleLights.size()  ) {

            //cout<<cpmap.storedPhotonNum<<endl;

			float offset_x = mRandom();
			float offset_y = mRandom();

			vec3 pos = simpleLights[i]->pos  - simpleLights[i]->size1*simpleLights[i]->dir1 *0.5+  simpleLights[i]->size1*simpleLights[i]->dir1*offset_x
								- simpleLights[i]->size2*simpleLights[i]->dir2*0.5   +  simpleLights[i]->size2 * simpleLights[i]->dir2*offset_y;
			//vec3 pos = simpleLights[i]->pos;

			vec3 nl(0, -1.0 ,0);

			float r1= 2*PI*mRandom(), r2=mRandom(), r2s=sqrt(r2);

			vec3 w=nl, u= vec3(1.0,0,0)^w ;
			u = u.normalize();
			vec3 v =w^u;

			vec3 d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2));
			d = d.normalize();

			//cout<<"pos "<<pos[0]<<"  "<<pos[1]<<"  "<<pos[2]<<endl;
			//cout<<"dir "<<d[0]<<"  "<<d[1]<<"  "<<d[2]<<endl;

			vec3 color(0.2,0.2,0.2);

			PhotonTraceCaustics( cpmap, simpleMeshes , pos , d , 0 , 2.40 , color , false ,causticMeshes);

		}

         cout<<"finish generating caustic photons"<<endl;

    }

    pmap.balance();
    pmap.ScalePhotonPower(1.0f/((float)pmap.storedPhotonNum>=1.0?(float)pmap.storedPhotonNum:1.0) ,0,pmap.storedPhotonNum);

    cpmap.balance();
    cpmap.ScalePhotonPower(1.0f/((float)cpmap.storedPhotonNum>=1.0?(float)cpmap.storedPhotonNum:1.0) ,0,cpmap.storedPhotonNum);


    ColorImage image;
	int x, y;
	Pixel p={0,0,0};

	image.init(ceil(cam.resX),ceil(cam.resY));
	//tempImage.init(cam.resX,cam.resY);


    vec2 vv;
    int q1 = 0 ,q2=0;
    int counter=0;



	for(x =0;x<cam.resX;x++){

        cout<<"row "<<x<<endl;

        for(y=0;y<cam.resY;y++){

            //cout<<"col "<<y<<endl;
            int r,g,b;
            r=0;
            g=0;
            b=0;

            vec3 color(0,0,0);
            for (int aa = 0; aa < KERNEL_SIZE; aa++)
			{
                vv[0] = x + AAFilter[aa][0];
                vv[1] = y + AAFilter[aa][1];

                 vec3 pixelCenter = cam.GetPixelCenter(vv);
                 vec3 ray = (pixelCenter-cam.center).normalize();

                 if(!TestVisualizeLight(cam.center,ray)){
                    color = color + AAFilter[aa][2] * RayCastSynthesisColorImage(cam,vv , pmap, cpmap );
                 }
                 else{
                    color = vec3(1.0,1.0,1.0);
                    rayHitAnything = true;
                    break;
                 }
            }

            //vec3 color = RayCastSynthesisColorImage(cam,vv , pmap, cpmap );
            //cout<<"color "<<color[0]<<" "<<color[1]<<" "<<color[2]<<endl;
            //averageCheck += (totalIntersectionCheck/(double)checkTime) / (double)(cam.resX*cam.resY);

            if(rayHitAnything){
                //cout<<"yooo"<<endl;
                //cout<<"draw white"<<q1++<<endl;
                //cout<<"hitanything  "<<counter++<<"  color:"<<color[0]<<" "<<color[1]<<" "<<color[2]<<endl;
                p.R = color[0] > 1.0f? 255 : (color[0]*255.0f);
                p.G = color[1] > 1.0f? 255 : (color[1]*255.0f);
                p.B = color[2] > 1.0f? 255 : (color[2]*255.0f);
                //cout<<"color "<<(int)p.R<<" "<<(int)p.G<<" "<<(int)p.B<<endl;
            }
            else{
                //cout<<"draw black"<<q2++<<endl;
                p.R = backGroundColor[0];
                p.G = backGroundColor[1];
                p.B = backGroundColor[2];
            }
            rayHitAnything = false;

            image.writePixel(x,y,p);

        }

	}

	//cout<<"yoyoy"<<endl;
	char s[100] = "colorOutput.ppm";

	image.outputPPM(s);

}

vec3 SimpleMeshRayTrace( vec3 ori , vec3 dir , int depth ,Camera cam , float currentNr){
    //cout<<"fuck"<<endl;
    return vec3(0,0,0);
    /* I don't implement normal ray trace in this build */
}

vec3 RayCastSynthesisColorImage(Camera  cam , vec2 xy , PhotonMap& pmap, PhotonMap& cpmap){
    vec3 pixelCenter = cam.GetPixelCenter(xy);
    vec3 ray = (pixelCenter-cam.center).normalize();

    return SimpleMeshRayTrace(pmap , cpmap, cam.center , ray , 0 , cam, 2.4 );
}


vec3 RayCastSynthesisColorImage(Camera  cam , vec2 xy){
    cout<<"fuuuck"<<endl;
    vec3 pixelCenter = cam.GetPixelCenter(xy);
    vec3 ray = (pixelCenter-cam.center).normalize();

    return SimpleMeshRayTrace(cam.center , ray , 0 , cam, airReractivity );
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


void GenerateColorImage(Camera cam , char* inputName){
    ColorImage image,tempImage;
	int x, y;
	Pixel p={0,0,0};

	image.init(ceil(cam.resX/cam.ssaaRatio),ceil(cam.resY/cam.ssaaRatio));
	tempImage.init(cam.resX,cam.resY);


    vec2 vv;
    int q1 = 0 ,q2=0;
    int counter=0;

     clock_t beg = clock();

	for(x =0;x<cam.resX;x++){
        cout<<"now render row:"<<x<<endl;
        for(y=0;y<cam.resY;y++){

            totalIntersectionCheck = 0;
            checkTime = 0;

            vv[0] = x;
            vv[1] = y;
            vec3 color = RayCastSynthesisColorImage(cam,vv);
            averageCheck += (totalIntersectionCheck/(double)checkTime) / (double)(cam.resX*cam.resY);

            if(rayHitAnything){
                //cout<<"draw white"<<q1++<<endl;
                //cout<<"hitanything  "<<counter++<<"  color:"<<color[0]<<" "<<color[1]<<" "<<color[2]<<endl;
                p.R = color[0] > 1.0f? 255 : (color[0]*255.0f);
                p.G = color[1] > 1.0f? 255 : (color[1]*255.0f);;
                p.B = color[2] > 1.0f? 255 : (color[2]*255.0f);;
            }
            else{
                //cout<<"draw black"<<q2++<<endl;
                p.R = backGroundColor[0];
                p.G = backGroundColor[1];
                p.B = backGroundColor[2];
            }
             rayHitAnything = false;
            tempImage.writePixel(x,y,p);
        }
	}

      clock_t ed = clock();
      elapsedTime = double(ed - beg) / CLOCKS_PER_SEC;
      printf("elapse time:%lf\n", elapsedTime);


	for(x =0;x<cam.resX/cam.ssaaRatio;x++){
        for(y=0;y<cam.resY/cam.ssaaRatio;y++){
            int r,g,b;
            r=0;
            g=0;
            b=0;
            //cout<<"x "<<x<<" y "<<y<<endl;
            for(int i=0;i<cam.ssaaRatio;i++){
                for(int j=0;j<cam.ssaaRatio;j++){
                    int temp1 = (x*cam.ssaaRatio +i) >= cam.resX ? cam.resX-1 :x*cam.ssaaRatio +i;
                    int temp2 = (y*cam.ssaaRatio +j) >= cam.resY ? cam.resY-1 :y*cam.ssaaRatio +j;

                    r += (int)tempImage.readPixel(temp1,temp2).R;
                    g += (int)tempImage.readPixel(temp1,temp2).G;
                    b += (int)tempImage.readPixel(temp1,temp2).B;
                }
            }
            //cout<<"r"<<r<<" g"<<g<<" b"<<b<<endl;
            r = r / (cam.ssaaRatio*cam.ssaaRatio);
            g = g / (cam.ssaaRatio*cam.ssaaRatio);
            b = b / (cam.ssaaRatio*cam.ssaaRatio);
            p.R = r;
            p.G = g;
            p.B = b;

            image.writePixel(x,y,p);

        }
	}

	//cout<<"yoyoy"<<endl;
	char s[100] = "colorOutput.ppm";

	image.outputPPM(s);

}


void GenerateBinaryIntersectionTestImage(Camera cam){

	ColorImage image;
	int x, y;
	Pixel p={0,0,0};

	image.init(cam.resX, cam.resY);
	//tempImage.init(cam.resX*cam.ssaaRatio , cam.resY * cam.ssaaRatio);

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

void MeasureRenderTime(Camera cam , char* inputName){

  //clock_t beg = clock();



  GenerateColorImage(cam,inputName);



  //clock_t ed = clock();
  //double elapsed_secs = double(ed - beg) / CLOCKS_PER_SEC;
  //printf("elapse time:%lf\n", elapsed_secs);

  char temp[100] = "executionTime_";

  ofstream fout(strcat(temp,inputName));

  if(!fout) {
        cout << "file output failed!" << endl;
        return;
  }

  fout<<"resolution:"<<cam.resX<<"X"<<cam.resY<<endl;
  fout<<"SSAA ratio:"<<cam.ssaaRatio<<endl;
  fout<<"elapse time:"<<elapsedTime<<" seconds"<<endl;
  fout<<"total kd-tree node count:"<<builded<<endl;
  fout<<"max depth of kd-tree subdivision:"<<maxDepth<<endl;
  fout<<"average ray-mesh Intersection check count:"<<averageCheck<<endl;
  fout.close();



}
