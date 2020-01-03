#include "GMath.h"
#include "GPoint.h"
#include "GPaint.h"
#include "GPixel.h"
#include "GMatrix.h"
#include <cmath>
		//  SX, KX, TX,
        //  KY, SY, TY,
		//fMat[0] = a;    fMat[1] = b;    fMat[2] = c;
       // fMat[3] = d;    fMat[4] = e;    fMat[5] = f;
        void GMatrix::setIdentity() {
			//GMatrix identity;
			this->set6(1, 0, 0, 0, 1, 0);
			//matrixStack.push_back(identity);
		}
        void GMatrix::setTranslate(float tx, float ty) {
			//GMatrix translate;
			this->set6(1, 0, tx, 0, 1, ty);
			//matrixStack.push_back(translate); 
		}
       void GMatrix::setScale(float sx, float sy) {
			//GMatrix scale;
			this->set6(sx, 0, 0, 0, sy, 0);
			//matrixStack.push_back(scale);
		}
       void GMatrix::setRotate(float radians){
			//GMatrix rotate;
			this->set6(cos(radians), -sin(radians), 0, sin(radians), cos(radians), 0);
			//matrixStack.push_back(rotate);
		}
        void GMatrix::setConcat(const GMatrix& secundo, const GMatrix& primo) {
			//Pts' = Secundo * Primo * Pts
			//secundo * primo -> [0][0]+ [1][3], [0][1] + [1][4], [0][2] +[1][5] +secu[2],
			// [3][0] +[4][3], [3][1] + [4][4], [3][2] + [4][5] + secu[5]
			 //Aa+Bd  Ab+Be  Ac+Bf+C 
			 //Da+Ed  Db+Ee  Dc+Ef+F
			//GMatrix concat;
			this->set6((secundo.fMat[0]*primo.fMat[0]) + (secundo.fMat[1]*primo.fMat[3]), ((secundo.fMat[0]*primo.fMat[1]) + (secundo.fMat[1]*primo.fMat[4])), (((secundo.fMat[0]*primo.fMat[2]) + (secundo.fMat[1]*primo.fMat[5])) + secundo.fMat[2]),
								((secundo.fMat[3]*primo.fMat[0]) + (secundo.fMat[4]*primo.fMat[3])), (secundo.fMat[3]*primo.fMat[1])+(secundo.fMat[4]*primo.fMat[4]), ((secundo.fMat[3]*primo.fMat[2]) + (secundo.fMat[4]*primo.fMat[5]) + secundo.fMat[5]));
			//matrixStack.push_back(concat);
		}
        bool GMatrix::invert(GMatrix* inverse) const {
			//We define the determinant of M to be: det = AE-BD
			//E -B  BF-EC
			//-D  A  DC-AF
			//*1/det
			float det = (this->fMat[0]*this->fMat[4])-(this->fMat[1]*this->fMat[3]);
			if(det == 0){
				return false;
			}
			//GMatrix inv;
			inverse->set6(this->fMat[4]/det, -this->fMat[1]/det, ((this->fMat[1]*this->fMat[5])-(this->fMat[4]*this->fMat[2]))/det, 
							-this->fMat[3]/det, this->fMat[0]/det, ((this->fMat[3]*this->fMat[2])-(this->fMat[0]*this->fMat[5]))/det);
			//matrixStack.push_back(inv);
			return true;				 				
		}
        void GMatrix::mapPoints(GPoint dst[], const GPoint src[], int count) const {
			//this.matrix * pts -> x = Ax + By + C, y = Dx+Ey+F
			/**
     *  Transform the set of points in src, storing the resulting points in dst, by applying this
     *  matrix. It is the caller's responsibility to allocate dst to be at least as large as src.
     *
     *  Note: It is legal for src and dst to point to the same memory (however, they may not
     *  partially overlap). Thus the following is supported.
     *
     *  GPoint pts[] = { ... };
     *  matrix.mapPoints(pts, pts, count);
     */
			float x;
			float y;
			GPoint dest;
			GPoint point;
			for(int i = 0; i < count; i++){
				point = src[i];
				x = (this->fMat[0]*point.x()) + (this->fMat[1]*point.y()) + this->fMat[2];
				y = (this->fMat[3]*point.x()) + (this->fMat[4]*point.y()) + this->fMat[5];
				dest.set(x, y);
				dst[i] = dest;
			}
		}