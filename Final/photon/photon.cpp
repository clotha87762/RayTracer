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


#define swap(ph,a,b) {Photon *ph2 = ph[a]; ph[a] =ph[b]; ph[b]=ph2;}

using namespace std;

float cosThetaTable[256];
float sinThetaTable[256];
float cosPhiTable[256];
float sinPhiTable[256];

vec3 GetPhotonDirection(Photon* p){

    vec3 temp;
    temp[0] = sinThetaTable[p->theta]*cosPhiTable[p->phi];
	temp[1] = sinThetaTable[p->theta]*sinPhiTable[p->phi];
    temp[2] = cosThetaTable[p->theta];
    return temp;
}



void InitDirectionTable(){

    for (int i = 0; i<256; i++) {
		double angle = double(i)*(1.0/256.0)*PI;
		cosThetaTable[i] = cos(angle);
		sinThetaTable[i] = sin(angle);
		cosPhiTable[i] = cos(2.0*angle);
		sinPhiTable[i] = sin(2.0*angle);
    }

}


 void PhotonMap::LocatePhoton(NearestPhotons *const np,const int index) const{

    Photon *p = &photons[index];
	float dist1;

	if (index<half_stored_photons) {
		dist1 = np->pos[p->plane] - p->position[p->plane];

		if (dist1>0.0) {
			LocatePhoton(np,2*index+1);
			if (dist1*dist1<np->dist2[0])
				LocatePhoton(np,2*index);
		} else {
			LocatePhoton(np,2*index);
			if (dist1*dist1 <np->dist2[0])
				LocatePhoton(np,2*index+1);
		}
	}

	dist1 = p->position[0] - np->pos[0];
	float dist2 = dist1*dist1;
	dist1 = p->position[1] - np->pos[1];
	dist2 += dist1*dist1;
	dist1 = p->position[2] - np->pos[2];
	dist2 += dist1*dist1;

	if (dist2<np->dist2[0]) {
		// we found a photon, Insert it in the candidate list
		if (np->found< np->_max) {
			// heap is not full; use array
			np->found++;
			np->dist2[np->found] = dist2;
			np->index[np->found] = p;
		} else {
			int j,parent;
			if (np->got_heap==0) { // need to build the heap
				float dst2;
                Photon *phot;
				int half_found = np->found>>1;
				for (int k=half_found; k>=1;k--){
					parent = k;
					phot = np->index[k];
					dst2 = np->dist2[k];
					while (parent <= half_found) {
						j = parent + parent;
						if (j<np->found && np->dist2[j]<np->dist2[j+1])
							j++;
						if (dst2>=np->dist2[j])
							break;
						np->dist2[parent] = np->dist2[j];
						np->index[parent] = np->index[j];
						parent = j;
					}
					np->dist2[parent] = dst2;
					np->index[parent] = phot;
				}
				np->got_heap = 1;
			}
			parent = 1;
			j = 2;
			while (j<=np->found) {
				if (j<np->found && np->dist2[j]<np->dist2[j+1])
					j++;
				if (dist2 > np->dist2[j])
					break;
				np->dist2[parent] = np->dist2[j];
				np->index[parent] = np->index[j];
				parent = j;
				j+=j;
			}
			np->index[parent] = p;
			np->dist2[parent] = dist2;

			np->dist2[0] = np->dist2[1];
		}
    }


}

void PhotonMap::balance(void){

    if (storedPhotonNum>1) {
		Photon **pa1 = (Photon**) malloc(sizeof(Photon*)*(storedPhotonNum+1));
		Photon **pa2 = (Photon**) malloc(sizeof(Photon*)*(storedPhotonNum+1));

		for (int i = 0; i<=storedPhotonNum; i++)
			pa2[i] = &photons[i];

		balance_segment(pa1, pa2, 1, 1, storedPhotonNum);
		free(pa2);

		int d, j = 1, foo = 1;
		Photon foo_photon = photons[j];

		for (int i=1; i<=storedPhotonNum; i++) {
			d = pa1[j]-photons;
			pa1[j] = 0;
			if (d != foo)
				photons[j] = photons[d];
			else {
				photons[j] = foo_photon;
				if (i<storedPhotonNum) {
					for (;foo<=storedPhotonNum;foo++)
						if (pa1[foo] != 0)
							break;
					foo_photon = photons[foo];
					j = foo;
				}
				continue;
			}
			j = d;
		}
		free(pa1);
	}
    half_stored_photons = storedPhotonNum/2 - 1;
}

void PhotonMap::balance_segment(
            Photon **pbal,
            Photon **porg,
            const int index,
            const int start,
            const int ending)
{
    int median = 1;
	while ((4*median)<=(ending-start+1))
		median += median;
	if ((3*median) <= (ending-start+1)) {
		median += median;
		median += start - 1;
	} else
		median = ending-median +1;


	int axis = 2;
	if ((bbMax[0]-bbMin[0])>(bbMax[1]-bbMin[1]) &&
		(bbMax[0]-bbMin[0])>(bbMax[2]-bbMin[2]))
		axis = 0;
	else if ((bbMax[1]-bbMin[1])>(bbMax[2]-bbMin[2]))
		axis = 1;

	median_split (porg, start, ending, median, axis);
	pbal[ index ] = porg[median];
	pbal[ index ]->plane = axis;


	if (median>start) {
		if (start<median-1) {
			const float tmp = bbMax[axis];
			bbMax[axis] = pbal[index]->position[axis];
			balance_segment(pbal, porg, 2*index, start, median-1);
			bbMax[axis] = tmp;
		} else {
			pbal[2*index] = porg[start];
		}
	}
	if (median<ending) {
		if (median+1<ending) {
			const float tmp = bbMin[axis];
			bbMin[axis] = pbal[index]->position[axis];
			balance_segment(pbal,porg, 2*index+1, median+1, ending);
			bbMin[axis] = tmp;
		} else {
			pbal[2*index+1] = porg[ending];
		}
    }


}

void PhotonMap::median_split(
            Photon **p,
            const int start,
            const int ending,
            const int median,
            const int axis )
{
    int left = start;
	int right = ending;

	while (right>left) {
		const float v = p[right]->position[axis];
		int i=left-1;
		int j=right;
		for(;;){
			while (p[++i]->position[axis]<v)
				;
			while (p[--j]->position[axis]>v && j>left)
				;
			if (i>=j)
				break;
			swap(p,i,j);
		}
		swap(p,i,right);
		if (i>=median)
			right = i-1;
		if (i<=median)
			left=i+1;
    }

}


PhotonMap::PhotonMap(int maxPhotons){

    this->maxPhotons = maxPhotons;
    photons = (Photon*) malloc(sizeof(Photon)*(maxPhotons+3));
    storedPhotonNum = 0;

    bbMin = vec3(-1e9,-1e9,-1e9);
    bbMax = vec3(1e9,1e9,1e9);

    if(photons == NULL){
        cout<<"memory is not enough for photons"<<endl;
        exit(-1);
    }

}

PhotonMap::PhotonMap(){
    this->maxPhotons = MAX_PHOTON;
     photons = (Photon*) malloc(sizeof(Photon)*(MAX_PHOTON+3));
     storedPhotonNum = 0;

     bbMin = vec3(-1e9,-1e9,-1e9);
     bbMax = vec3(1e9,1e9,1e9);

      if(photons == NULL){
        cout<<"memory is not enough for photons"<<endl;
        exit(-1);
    }
}

PhotonMap::~PhotonMap(){

    free(photons);
    cout<<"photon map is deallocated"<<endl;

}

SimplePhotonQueryResult PhotonMap::SimpleQueryPhotons(const SimplePhotonQuery& query){

    SimplePhotonQueryResult result;
    vec3 pos = query.pos;
    int numPhoton = query.numPhoton;
    float maxDis = query.maxDist;

    /* 一定要記得deallocate啊啊啊*/

    //float* disArray = (float*) malloc(sizeof(float)*numPhoton);
    //int* indexArray = (int*) malloc(sizeof(int) * numPhoton );

    //SimplePhotonQueryResult result;

    result.dist = (float*) malloc(sizeof(float)*numPhoton);
    result.indexes = (int*) malloc(sizeof(int) * numPhoton );
    result.found = 0;

    vector<IndexDistPair> pairs;

    /* 極慢....*/
    for(int i=0;i<storedPhotonNum;i++){

        const Photon* node = &photons[i];


        float dis = vec3(node->position - pos).length(); // 只是比較 避免浪費時間算sqrt

        if(dis > maxDis)
            continue;

        IndexDistPair newPair;
        newPair.dist = dis;
        newPair.index = i;
        pairs.push_back(newPair);
    }

    std::sort( pairs.begin(),pairs.end());

    for(int i=0;i<numPhoton && i < pairs.size();i++){
        result.dist[i] = sqrt( pairs[i].dist );
        result.indexes[i] = pairs[i].index;
        result.found++;
    }

    return result;
}

vec3 PhotonMap::CalcIrradiance(vec3 pos, vec3 normal,  float maxDistance ,  int photonNum ,  int minPhotons ){

    vec3 irradiance(0,0,0);

    if(QUERY_MODE == NAIVE){

        SimplePhotonQuery query;
        query.numPhoton = photonNum;
        query.maxDist = maxDistance;
        query.pos = pos;

        float radius = 1.0f;

        SimplePhotonQueryResult result = SimpleQueryPhotons(query);

        if(result.found < minPhotons)
            return vec3(0,0,0);

        //vec3 irradiance(0,0,0);

        for (int i=0;i<result.found; i++){

            Photon *p = &photons[result.indexes[i]];

            vec3 pdir = GetPhotonDirection(p);

            if ((pdir[0]*normal[0]+pdir[1]*normal[1]+pdir[2]*normal[2])<0.0f) {

                for(int j=0;j<3;j++)
                    irradiance[j] = irradiance[j] + p->power[j];
            }

        }

        if(result.found>0){
            radius = result.dist[result.found-1];
        }

        irradiance *=  1.0 / ( radius*radius*PI );

        //cout<<"irradiance "<<irradiance[0]<<" "<<irradiance[1]<<" "<<irradiance[2]<<endl;

        // avoid memory leak
        free(result.dist);
        free(result.indexes);

        return irradiance;

    }
    else if(QUERY_MODE == BALANCE_TREE){

        //cout<<"whattt"<<endl;
        vec3 irrad(0,0,0);
        //irrad[0] = irrad[1] = irrad[2] = 0.0;

        NearestPhotons np;
        np.dist2 = (float*)malloc(sizeof(float)*(photonNum+1));
        np.index = (Photon**)malloc(sizeof(Photon*)*(photonNum+1));

        np.pos[0] = pos[0]; np.pos[1] = pos[1]; np.pos[2] = pos[2];
        np._max = photonNum;
        np.found = 0;
        np.got_heap = 0;
        np.dist2[0] = maxDistance*maxDistance;

        LocatePhoton(&np,1);

        //cout<<"found"<<np.found<<endl;

        if (np.found<8){
            free(np.dist2);
            free(np.index);
            return vec3(0,0,0);
        }
        //float pdir[3];
        vec3 pdir;
        for (int i=1;i<=np.found; i++){

            Photon *p = np.index[i];
            pdir = GetPhotonDirection(p);

            if ((pdir[0]*normal[0]+pdir[1]*normal[1]+pdir[2]*normal[2])<0.0f) {
                irrad[0] += p->power[0];
                irrad[1] += p->power[1];
                irrad[2] += p->power[2];
            }

        }

        const float tmp = (1.0f/PI)/(np.dist2[0]);

        irrad *= tmp;

        //cout<<"irradiance "<<irrad[0]<<"  "<<irrad[1]<<"  "<<irrad[2]<<endl;

        free(np.dist2);
        free(np.index);

        return irrad;
    }

}

void PhotonTrace( PhotonMap& photonMap,vector<SimpleMesh*>& simpleMeshes, vec3 ori, vec3 dir, int depth , float currentNr , vec3 lightColor){

    if(depth>=PHOTON_TRACE_MAX_DEPTH)
        return ;

    float nearestDis = 1e9;

    SimpleMesh* nearestMesh = 0;
    vector<SimpleMesh*>::iterator meshIte;

    vec3 colour(0,0,0);
    vec3 hitPoint;

    //checkTime++;

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



    vec3 ambient(0,0,0);
    vec3 diffuse(0,0,0);
    vec3 specular(0,0,0);



    if(nearestMesh){

        //nearestMesh->PrintElements();

        vec3 normal = nearestMesh->GetNormal(hitPoint);
        //vec3
        vec3 normalL = vec3(normal);
        // 問題 該不會出在normal要不要轉正上吧 qq

        if(normal*dir > 0 ){
            normalL = normalL* -1.0f;
        }

        vec3 reflected;
        vec3 transmitted;

        float REFLECT_PROB = nearestMesh->reflect;
        float REFRACT_PROB = nearestMesh->refract;

        float ABSORP_PROB = 1.0 - REFLECT_PROB - REFRACT_PROB > 0? 1.0 - REFLECT_PROB - REFRACT_PROB : 0;

        float russianRoulette = mRandom();

        if(russianRoulette < REFLECT_PROB && depth < PHOTON_TRACE_MAX_DEPTH ){  // reflect photon with full energy
            float dot = (dir*normal);
            reflected = dir - 2 * dot * normal;
            float wd = nearestMesh->Kd;
            float ws = nearestMesh->Ks;

            russianRoulette = mRandom() * (wd+ws);

            // specularly reflect the photon
            PhotonTrace(  photonMap, simpleMeshes, hitPoint , reflected, depth +1 , currentNr , lightColor);

        }
        else if(russianRoulette < REFLECT_PROB + REFRACT_PROB  && depth < PHOTON_TRACE_MAX_DEPTH  ){ // refract photon with full energy

            vec3 reflRay( dir - normal*2*normal*dir);
            bool into = normal * normalL >0;
            float nc=1, nt=currentNr , nnt=into?nc/nt:nt/nc, ddn= dir * normalL , cos2t;
            if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    {

                PhotonTrace( photonMap,simpleMeshes , hitPoint , reflRay  , depth+1 , currentNr ,  lightColor );
            }

            vec3  tdir = (dir*nnt - normal *((into?1:-1)*(ddn*nnt+sqrt(cos2t))));
            tdir = tdir.normalize();

            float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn: tdir*normal);

            float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re, RP=Re/P, TP=Tr/(1-P);
            if(russianRoulette > REFLECT_PROB ){
                // Russian roulette
                PhotonTrace( photonMap ,simpleMeshes, hitPoint ,reflRay , depth+1 ,currentNr, lightColor * RP );

            }
            else{
                PhotonTrace (photonMap , simpleMeshes  ,hitPoint , tdir , depth+1 ,currentNr ,lightColor * TP );
            }
            /*
            float n1 = currentNr;
            float n2 = nearestMesh->Nr;
            float cos1 = dir * normal;
            float sin2 = (n1/n2) * sqrt( 1-(cos1*cos1));
            transmitted = (   ((n1/n2)*dir) + ( cos1*(n1/n2)  - sqrt(1-(sin2*sin2)))*normal   ).normalize();
            vec3 displacePoint = hitPoint - normal * refractBias;
            */

            //PhotonTrace(  photonMap,  displacePoint , transmitted , depth +1 , currentNr , lightColor);

        }
        else{  //Absorption
            vec3 sc = nearestMesh->surfColor;
            vec3 power ( sc[0]*lightColor[0] / ( ABSORP_PROB )  ,   sc[1]*lightColor[1] / (ABSORP_PROB )  , sc[2]*lightColor[2] / ( ABSORP_PROB ) );
            photonMap.AddPhoton(hitPoint , dir , power);
            return;
        }


    }


    return;

}

void PhotonTraceCaustics(PhotonMap& photonMap,vector<SimpleMesh*>& simpleMeshes , vec3 ori , vec3 dir, int depth, float currentNr , vec3 lightColor , bool flag ,vector<SimpleMesh*>& cMeshes){

    // 強制他打在地板上xd
    if(dir[1] > 0)
        return;

    if (!flag){
        bool bb =false;
        for(int i=0;i<cMeshes.size();i++){

         bool b = cMeshes[i]->IntersectionTest(dir,ori);
             if(b){
                bb = true;
                break;
             }
        }
        if(!bb)
            return;
    }

    //cout<<"ori "<<ori[0]<<"  "<<ori[1]<<"  "<<ori[2]<<endl;
    //cout<<"dir "<<dir[0]<<"  "<<dir[1]<<"  "<<dir[2]<<endl;


    if(depth>=2)
        return ;

    float nearestDis = 1e9;

    SimpleMesh* nearestMesh = 0;
    vector<SimpleMesh*>::iterator meshIte;

    vec3 colour(0,0,0);
    vec3 hitPoint;



    //checkTime++;

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



    vec3 ambient(0,0,0);
    vec3 diffuse(0,0,0);
    vec3 specular(0,0,0);



    if(nearestMesh){

        if (flag==false && nearestMesh->refract < 0.01) return;
        flag = true;

        //nearestMesh->PrintElements();

        vec3 normal = nearestMesh->GetNormal(hitPoint);
        //vec3
        vec3 normalL = vec3(normal);
        // 問題 該不會出在normal要不要轉正上吧 qq

        if(normal*dir > 0 ){
            normalL = normalL* -1.0f;
        }

        vec3 reflected;
        vec3 transmitted;

        float REFLECT_PROB = nearestMesh->reflect;
        float REFRACT_PROB = nearestMesh->refract;

        float ABSORP_PROB = 1.0 - REFLECT_PROB - REFRACT_PROB > 0? 1.0 - REFLECT_PROB - REFRACT_PROB : 0;

        float russianRoulette = mRandom();
        /*
        if(russianRoulette < REFLECT_PROB && depth < PHOTON_TRACE_MAX_DEPTH ){  // reflect photon with full energy
            float dot = (dir*normal);
            reflected = dir - 2 * dot * normal;
            float wd = nearestMesh->Kd;
            float ws = nearestMesh->Ks;

            russianRoulette = mRandom() * (wd+ws);

            // specularly reflect the photon
            PhotonTraceCaustics(  photonMap, simpleMeshes, hitPoint , reflected, depth +1 , currentNr , lightColor  ,flag ,cMeshes);

        }
        else*/
        if(russianRoulette <  REFRACT_PROB  && depth < PHOTON_TRACE_MAX_DEPTH  ){ // refract photon with full energy

            vec3 reflRay( dir - normal*2*normal*dir);
            bool into = normal * normalL >0;
            float nc=1, nt=currentNr , nnt=into?nc/nt:nt/nc, ddn= dir * normalL , cos2t;
            if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    {

                //PhotonTraceCaustics( photonMap,simpleMeshes , hitPoint , reflRay  , depth+1 , currentNr ,  lightColor ,flag ,cMeshes);
            }

            vec3  tdir = (dir*nnt - normal *((into?1:-1)*(ddn*nnt+sqrt(cos2t))));
            tdir = tdir.normalize();

            float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn: tdir*normal);

            float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re, RP=Re/P, TP=Tr/(1-P);
            PhotonTraceCaustics (photonMap , simpleMeshes  ,hitPoint , tdir , depth+1 ,currentNr ,lightColor  ,flag ,cMeshes);

            //PhotonTrace(  photonMap,  displacePoint , transmitted , depth +1 , currentNr , lightColor);

        }
        else if(flag){  //Absorption
            vec3 sc = nearestMesh->surfColor;
            vec3 power ( sc[0]*lightColor[0]   ,   sc[1]*lightColor[1]   , sc[2]*lightColor[2] );

            if(hitPoint[1]>-0.3)
              return;

            photonMap.AddPhoton(hitPoint , dir , power);
            return;
        }


    }


    return;

}



vec3 TraceIrradiance(PhotonMap& photonMap, vector<SimpleMesh*>& simpleMeshes, vec3 ori , vec3 dir , int depth , float currentNr){



    if(depth>IRRADIANCE_MAX_DEPTH)
        return vec3(0,0,0);

    float nearestDis = 1e9;

    SimpleMesh* nearestMesh = 0;
    vector<SimpleMesh*>::iterator meshIte;

    vec3 colour(0,0,0);
    vec3 hitPoint;

    //checkTime++;

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

        bool isGoOut = false;


        //nearestMesh->PrintElements();

        vec3 normal = nearestMesh->GetNormal(hitPoint);
        vec3 normalL = normal;
        if(normal*dir > 0 ){
            normalL = normal* -1.0f;
            isGoOut = true;
        }


        vec3 reflected;
        vec3 transmitted;
        //vec3 col = photonMap.CalcIrradiance(hitPoint, normal ,  QUERY_MAX_DIS,  PHOTON_QUERY_COUNT , MIN_PHOTON_TO_ESTIMATE );
        //colour = colour + col;



        if(nearestMesh->Kd > 0){

            vec3 col =  photonMap.CalcIrradiance(hitPoint, normal ,  QUERY_MAX_DIS,  PHOTON_QUERY_COUNT , MIN_PHOTON_TO_ESTIMATE );

            colour = colour +  col * nearestMesh->Kd;
        }

        if(nearestMesh->reflect > 0.01f){

            float dot = (dir*normal);
            reflected = dir - 2 * dot * normal;

            vec3 reflectColor =  TraceIrradiance(photonMap ,simpleMeshes, hitPoint , reflected ,  depth+1 , currentNr);//SimpleMeshRayTrace(hitPoint,reflected ,depth+1,cam ,currentNr ); // reflected是否要反過來???
            colour = colour + reflectColor * nearestMesh->reflect ;
        }

        if(nearestMesh->refract > 0.01f){

            vec3 reflRay(dir- ((normal*2*normal) * dir));
            bool into = normal * normalL >0;
            float nc=1, nt=currentNr , nnt = into?nc/nt:nt/nc, ddn=dir*normalL, cos2t;
            if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    {// Total internal reflection
                //return  TraceIrradiance(photonMap,simpleMeshes,hitPoint,reflRay,depth+1 , currentNr);
            }

            vec3 tdir = (dir*nnt - normal*((into?1:-1)*(ddn*nnt+sqrt(cos2t))));
            tdir = tdir.normalize();
            float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn: tdir*normal );
            float Re=R0+(1-R0)*c*c*c*c*c;
            float Tr=1-Re;
            return nearestMesh->emissionColor + nearestMesh->refract * (  TraceIrradiance(photonMap,simpleMeshes , hitPoint , reflRay ,  depth+1 , currentNr) *Re  +
                                                                          TraceIrradiance(photonMap,simpleMeshes ,hitPoint , tdir ,depth+1, currentNr)*Tr);

        }



    }


    return colour;

}





void PhotonMap::ScalePhotonPower(float scale,int start,int num){

    for(int i=start;i<start+num ;i++){
        for(int j=0;j<3;j++)
            photons[i].power[j] =  photons[i].power[j] * scale;
    }

}

void PhotonMap::AddPhoton(vec3 pos, vec3 inDir, vec3 power){

    if(storedPhotonNum + 1 > maxPhotons)
        return;


    //強制他打在地板上xd
    //if(pos[1]>-0.3)
    //   return;

    //pos = pos.normalize();
    inDir = inDir.normalize();
    //power = power.normalize();
    //cout<<"pos "<<pos[0]<<"  "<<pos[1]<<"  "<<pos[2]<<endl;



    Photon* toStored = &photons[++storedPhotonNum];

    for(int i=0;i<3;i++){

        if(pos[i] > bbMax[i])
            bbMax[i] = pos[i];
        if(pos[i] < bbMin[i])
            bbMin[i] = pos[i];

    }

    //cout<<"power "<<power[0]<<"  "<<power[1]<<"  "<<power[2]<<endl;

    toStored->position = vec3(pos[0],pos[1],pos[2]);
    toStored->power = vec3(power[0],power[1],power[2]);
    toStored->inDir = vec3(inDir[0],inDir[1],inDir[2]);

    int phi =   (int) ( atan2(inDir[1],inDir[0])*(256.0/(2.0*M_PI)) );
    int theta = (int)  (acos(inDir[2]) * (256.0 / PI));

    if(phi > 255)
        phi = 255;

    if(phi < 0)
        phi = (unsigned char) (phi + 256);

    toStored->phi = (unsigned char) phi;

    if(theta > 255)
        theta = 255;
    toStored->theta = (unsigned char) theta;


}
