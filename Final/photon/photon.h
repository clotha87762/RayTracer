#ifndef __photon_h__
#define __photon_h__

#include <random>
#include <ctime>
#include <functional>

#include "algebra3.h"
#include "parameter.h"
#include "Mesh.h"
#include "Camera.h"
//#include "final.h"

using namespace std;

static default_random_engine generator(time(nullptr));
static uniform_real_distribution<float> distribution(0,1);
static auto mRandom = std::bind(distribution,generator);


/*
void InitRandom(){
     srand( time(NULL) );
}

float mRandom(){
    return (float) rand() / (RAND_MAX + 1.0);
}
*/

class SimpleLight{

    public :

    float intensity;
    vec3 pos;
    vec3 color;

    vec3 dir1, dir2;
    int sampleCountDir1, sampleCountDir2;
    float size1,size2;

    SimpleLight() : intensity(1.0f) , pos(0,0,0) , color(1.0f,1.0f,1.0f){

    }

    virtual vector<vec3> GetSamplePoses(){
        vector<vec3> temp;
        temp.push_back(pos);
        return temp;
    };



};



class SpotLight : public SimpleLight{
    public :
    vec3 lightDir;
    float angle;
    SpotLight() : SimpleLight(),lightDir(0,0,1.0f),angle(60.0f){
    }
};

class DirectionalLight : public SimpleLight{
    public :

    vec3 lightDir;
    DirectionalLight() : SimpleLight(),lightDir(0,0,1.0f){

    }


};


class RectLight : public SimpleLight{

    public:

    //vec3 dir1, dir2;
   // int sampleCountDir1, sampleCountDir2;
    //float size1,size2;

    RectLight() : SimpleLight(){
    }

    vector<vec3> GetSamplePoses(){

        dir1 = dir1.normalize();
        dir2 = dir2.normalize();

        vec3 startPos = pos - (dir1*size1/2.0) - (dir2*size2/2.0);
        vec3 endPos = pos + (dir1*size1/2.0) + (dir2*size2/2.0);

        vector<vec3> samples;
        //vec3 start =
        vec3 nowPos(0,0,0);
        for(int i=0;i<sampleCountDir1;i++){
            nowPos = startPos + dir1 * size1 *i/ (float)sampleCountDir1;
            for(int j=0;j<sampleCountDir2;j++){
                vec3 tempPos = nowPos + dir2 * size2 * j / (float)sampleCountDir2;
                samples.push_back(tempPos);
            }

        }

        return samples;
    }
};

class PointLight : public SimpleLight{
    public :
    PointLight() : SimpleLight(){
    }

    vector<vec3> GetSamplePoses(){
        vector<vec3> poses;
        poses.push_back(vec3(pos));
        return poses;
    }
};

class SphereLight : public SimpleLight{
    vec3 pos;
    float radius;
    float sampleDensity; // samples per area
};


struct PhotonStruct{

    vec3 position;
    vec3 inDir;
    vec3 power;
    unsigned char phi , theta;
    short plane;

};


struct IndexDistPair{

    public:

    int index;
    float dist;

    bool operator < (const IndexDistPair& idp) const
    {
        return (dist < idp.dist);
    }
};



struct SimplePhotonQueryResultStruct{
    int* indexes;
    float* dist;
    int found;
};

struct SimplePhotonQueryStruct{
    float maxDist;
    int numPhoton;
    vec3 pos;
};

struct PhotonQueryResultStruct{

};

struct PhotonQueryStruct{

};

typedef struct PhotonStruct Photon;
typedef struct PhotonQueryResultStruct PhotonQueryResult;
typedef struct PhotonQueryStruct PhotonQuery;
typedef struct SimplePhotonQueryResultStruct SimplePhotonQueryResult;
typedef struct SimplePhotonQueryStruct SimplePhotonQuery;




extern float cosThetaTable[256];
extern float sinThetaTable[256];
extern float cosPhiTable[256];
extern float sinPhiTable[256];
void InitDirectionTable();


typedef struct NearestPhotons {
	int _max;
	int found;
	int got_heap;
	float pos[3];
	float *dist2;
	Photon **index;
} NearestPhotons;


class PhotonMap{

    public:

        PhotonMap(int maxPhotons);
        PhotonMap();
        ~PhotonMap();

        void AddPhoton(vec3 pos, vec3 inDir, vec3 power);

        PhotonQueryResult QueryPhotons(const PhotonQuery& query);
        SimplePhotonQueryResult SimpleQueryPhotons(const SimplePhotonQuery& query);

        // return flux of RGB Channel
        vec3 CalcIrradiance(vec3 pos, vec3 normal,  float maxDistance ,  int photonNum ,  int minPhotons );
        void ScalePhotonPower(float scale,int start,int num);


        void balance(void);

        void LocatePhoton(
		NearestPhotons *const np,
        const int index) const;

        void balance_segment(
            Photon ** pbal,
            Photon ** porg,
            const int index,
            const int start,
            const int ending);

        void median_split(
            Photon ** p,
            const int start,
            const int ending,
            const int median,
            const int axis );


        int half_stored_photons;

        int prev_scale;

        Photon* photons;
        int maxPhotons;
        int storedPhotonNum;
        vec3 bbMax, bbMin;

};

//void InitDirectionTable();
void PhotonTrace(PhotonMap& photonMap,vector<SimpleMesh*>& simpleMeshes, vec3 ori, vec3 dir, int depth , float currentNr , vec3 lightColor);
void PhotonTraceCaustics(PhotonMap& photonMap, vector<SimpleMesh*>& simpleMeshes, vec3 ori , vec3 dir, int depth, float currentNr , vec3 lightColor , bool flag , vector<SimpleMesh*>& cmesh);
vec3 GetPhotonDirection( Photon* photon);

vec3 SimpleMeshRayTrace(PhotonMap& pmap,PhotonMap& cpmap, vec3 ori , vec3 dir , int depth ,Camera cam , float currentNr);
vec3 RayCastSynthesisColorImage(Camera cam, vec2 xy, PhotonMap& pamp, PhotonMap& cpmap);

vec3 TraceIrradiance(PhotonMap& photonMap ,vector<SimpleMesh*>& simpleMeshes, vec3 ori , vec3 dir , int depth , float currentNr);

#endif
