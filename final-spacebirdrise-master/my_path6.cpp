#include "GPath.h"
#include "GMatrix.h"
#include "GRect.h"
#include "GPoint.h"
#include <iostream>
#include <string>
#include <vector>
GPath& GPath::addRect(const GRect& r, Direction dir){
//moveto start a new contour
//lineto continue
GPoint p1;
p1.fX = r.left();
p1.fY = r.top();
GPoint p2;
p2.fX = r.right();
p2.fY = r.top();
GPoint p3;
p3.fX = r.right();
p3.fY = r.bottom();
GPoint p4;
p4.fX = r.left();
p4.fY = r.bottom();

if(dir == Direction::kCW_Direction){
	this->moveTo(p1).lineTo(p2).lineTo(p3).lineTo(p4);
}else{
	this->moveTo(p1).lineTo(p4).lineTo(p3).lineTo(p2);
}
//GPath path = *this;
//GPath& ref = path;
//GPath path;
//Verb verbs[4];
//verbs = Edger(path);
//for(int i = 0;  i < 4; i++){
//	cout << "KLine? " << verbs[i]; 
//}
GPoint pts[4];

return *this;
/*
    enum Direction {
        kCW_Direction,
        kCCW_Direction,
    };
	*/	
	//GBlendMode::kClear
}

GPath& GPath::addPolygon(const GPoint pts[], int count){
GPath path;
for(int i = 0; i < count; i ++){
	if(i == 0){
		this->moveTo(pts[i]);
	}else{
	this->lineTo(pts[i]);
	}
}
//path = *this;
//GPath& ref = path;
//GPath path;
return *this;
	
}
           
		  


GRect GPath::bounds() const{
GRect rect;
int left;
int right;
int top;
int bottom;
//GPath replace = this;

//this->reset();
//GPath its made up of edges, creating an edger, then edger.next(pts[]) will result in the pts array adding two new pts to the array in pts[0] and pts[1] and the action itself resolves 
//to a verb value
/*
enum Verb {
        kMove,  // returns pts[0] from Iter
        kLine,  // returns pts[0]..pts[1] from Iter
        kDone   // returns nothing in pts, Iter is done
    };
*/
GPoint p;
p.fX = 100000;
p.fY = 100000;
GPoint points[2] = {p};
Iter edger = Iter(*this);

Verb verbCheck = edger.next(points);
//std::cout << verbCheck;
if(points[0].x() == 100000 && points[0].y() == 100000){
	//std::cout << verbCheck;
	return GRect::MakeWH(0,0);
}
if(verbCheck == Verb::kDone){
	return GRect::MakeWH(0,0);
}else if(verbCheck == Verb::kMove){
	GPoint start = points[0];
	float x = points[0].x();
	float y = points[0].y();
	//std::cout << "X is " << x;
	//std::cout << "Y is " <<  y;
	left = start.x();
	top = start.y();
	right = start.x();
	bottom = start.y();
	rect.setXYWH(start.x(),start.y(),0,0);
	//return rect;
}

while(verbCheck != Verb::kDone){
	if(verbCheck == Verb::kLine){
	GPoint current = points[0];
	GPoint current2 = points[1];
	if(current.x() < left){
		left = current.x();
	}
	if(current.x() > right){
		right = current.x();
	}
	if(current.y() < top){
		top = current.y();
	}
	if(current.y() > bottom){
		bottom = current.y();
	}
	if(current2.x() < left){
		left = current2.x();
	}
	if(current2.x() > right){
		right = current2.x();
	}
	if(current2.y() < top){
		top = current2.y();
	}
	if(current2.y() > bottom){
		bottom = current2.y();
	}
	}
	verbCheck = edger.next(points);
}
rect.setLTRB(left, top, right, bottom);
//rect = rectFixer(rect)

return rect;
	
}




void GPath::transform(const GMatrix& m){
//GPath its made up of edges, creating an edger, then edger.next(pts[]) will result in the pts array adding two new pts to the array in pts[0] and pts[1] and the action itself resolves as a variable for a verb
for(int i = 0; i < this->fPts.size(); i++){
	GPoint newPoint;
	GPoint point = this->fPts[i];
	m.mapPoints(&newPoint, &point, 1);
	this->fPts[i] = newPoint;
}
return;
/* GPoint points[2];
 bool firstEdge = true;
 std::vector<GPoint> pointContainer;
 int index = 0;
 Edger edger = Edger(*this);
 Verb verbCheck = edger.next(points);
 while(verbCheck != Verb::kDone){
 
 if(verbCheck == Verb::kMove){
 	GPoint start = points[0];
	m.mapPoints(&start, &start, 1);
	pointContainer.push_back(start);
	index = index + 1;
 }else if(verbCheck == Verb::kLine){
 GPoint current = points[0];
 GPoint current2 = points[1];
 m.mapPoints(&current, &current, 1);
 m.mapPoints(&current2, &current2, 1);
 if(firstEdge){
 	pointContainer.push_back(current);
	pointContainer.push_back(current2);
	index = index + 2;
	firstEdge = false;
 }else{
 pointContainer.push_back(current2);
 index = index + 1;
 }
 }
 verbCheck = edger.next(points);
 }
 //this->reset();
 GPoint* pC = &pointContainer[0];
 //this->addPolygon(pC, index);
return;	
*/
}
GPath& GPath::addCircle(GPoint center, float radius, GPath::Direction dir){
GPoint point[16];
	//point[0] = GPoint::Make(0,0);
	point[0] = GPoint::Make(1,0);
	point[1] = GPoint::Make(1, tan(M_PI/8));
	point[2] = GPoint::Make(sqrt(2)/2, sqrt(2)/2);
	point[3] = GPoint::Make (tan(M_PI/8), 1);
	point[4] = GPoint::Make(0,1);
	point[5] = GPoint::Make(-tan(M_PI/8), 1);
	point[6] = GPoint::Make(-sqrt(2)/2, sqrt(2)/2);
	point[7] = GPoint::Make(-1, tan(M_PI/8));
	point[8] = GPoint::Make(-1, 0);
	point[9] = GPoint::Make(-1, -tan(M_PI/8));
	point[10] = GPoint::Make(-sqrt(2)/2, -sqrt(2)/2);
	point[11] = GPoint::Make(-tan(M_PI/8), -1);
	point[12] = GPoint::Make(0,-1);
	point[13] = GPoint::Make(tan(M_PI/8), -1);
	point[14] = GPoint::Make(sqrt(2)/2, -sqrt(2)/2);
	point[15] = GPoint::Make(1, -tan(M_PI/8));
	GMatrix m;
	GMatrix p1;
	GMatrix p2;
	p1.setTranslate(center.x(), center.y());
	p2.setScale(radius, radius);
	m.setConcat(p1,p2);
	m.mapPoints(point, 16);
	if(dir == GPath::Direction::kCW_Direction){
		this->moveTo(point[0]).quadTo(point[1], point[2]).quadTo(point[3], point[4]).quadTo(point[5], point[6]).quadTo(point[7],point[8]).quadTo(point[9],point[10]).quadTo(point[11],point[12]).quadTo(point[13],point[14]).quadTo(point[15],point[0]);
	}
	if(dir == GPath::Direction::kCCW_Direction){
		this->moveTo(point[0]).quadTo(point[15], point[14]).quadTo(point[13], point[12]).quadTo(point[11], point[10]).quadTo(point[9],point[8]).quadTo(point[7],point[6]).quadTo(point[5],point[4]).quadTo(point[3],point[2]).quadTo(point[1],point[0]);
	}
	
	return *this;
}
void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t){
	GPoint a = src[0];
	GPoint b = src[1];
	GPoint c = src[2];
	GPoint ab =  GPoint::Make(a.x()+(t*(b.x() - a.x())), a.y() + (t*(b.y() - a.y())));
	GPoint bc =  GPoint::Make(b.x() + (t*(c.x() - b.x())), b.y() + (t*(c.y() - b.y())));
	GPoint abc = GPoint::Make((a.x()+(t*(b.x() - a.x()))) + (t*((b.x() + (t*(c.x() - b.x()))) - (a.x()+(t*(b.x() - a.x()))))), (a.y()+(t*(b.y() - a.y()))) + (t*((b.y() + (t*(c.y() - b.y()))) - (a.y()+(t*(b.y() - a.y()))))));
	dst[0] = a;
	dst[1] = ab;
	dst[2] = abc;
	dst[3] = bc;
	dst[4] = c;
	return;
}
void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t){
	GPoint a = src[0];
	GPoint b = src[1];
	GPoint c = src[2];
	GPoint d = src[3];
	
	GPoint ab =  GPoint::Make(a.x()+(t*(b.x() - a.x())), a.y() + (t*(b.y() - a.y())));
	GPoint bc = GPoint::Make(b.x() + (t*(c.x() - b.x())), b.y() + (t*(c.y() - b.y())));
	GPoint cd = GPoint::Make(c.x() + (t*(d.x() - c.x())), c.y() + (t*(d.y() - c.y())));
	//abc = AB + t(BC - AB)
	GPoint abc = GPoint::Make((a.x()+(t*(b.x() - a.x()))) + (t*((b.x() + (t*(c.x() - b.x()))) - (a.x()+(t*(b.x() - a.x()))))), (a.y()+(t*(b.y() - a.y()))) + (t*((b.y() + (t*(c.y() - b.y()))) - (a.y()+(t*(b.y() - a.y()))))));
	//abcd = ABC + t(BCD - ABC)
	GPoint abcd = GPoint::Make(((a.x()+(t*(b.x() - a.x()))) + (t*((b.x() + (t*(c.x() - b.x()))) - (a.x()+(t*(b.x() - a.x())))))) + t*(((b.x() + (t*(c.x() - b.x()))) + (t*((c.x() + (t*(d.x() - c.x()))) - (b.x() + (t*(c.x() - b.x())))))) - ((a.x()+(t*(b.x() - a.x()))) + (t*((b.x() + (t*(c.x() - b.x()))) - (a.x()+(t*(b.x() - a.x()))))))), ((a.y()+(t*(b.y() - a.y()))) + (t*((b.y() + (t*(c.y() - b.y()))) - (a.y()+(t*(b.y() - a.y())))))) + t*(((b.y() + (t*(c.y() - b.y()))) + (t*((c.y() + (t*(d.y() - c.y()))) - (b.y() + (t*(c.y() - b.y())))))) - ((a.y()+(t*(b.y() - a.y()))) + (t*((b.y() + (t*(c.y() - b.y()))) - (a.y()+(t*(b.y() - a.y()))))))));
	//BCD = BC + t(CD - BC)
	GPoint bcd = GPoint::Make((b.x() + (t*(c.x() - b.x()))) + (t*((c.x() + (t*(d.x() - c.x()))) - (b.x() + (t*(c.x() - b.x()))))) ,(b.y() + (t*(c.y() - b.y()))) + (t*((c.y() + (t*(d.y() - c.y()))) - (b.y() + (t*(c.y() - b.y()))))) );
	
	dst[0] = a;
	dst[1] = ab;
	dst[2] = abc;
	dst[3] = abcd;
	dst[4] = bcd;
	dst[5] = cd;
	dst[6] = d;
}
