#include "GBitmap.h"
#include "GCanvas.h"
#include "GBlendMode.h"
#include "GColor.h"
#include "GMath.h"
#include "GPoint.h"
#include "GPaint.h"
#include "GPixel.h"
#include "GRect.h"
#include "GMatrix.h"
#include "GShader.h"
//#include "GFilter.h"
#include "GPath.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
struct Layer{
	GBitmap bounds;
	GMatrix ctm;
	GPaint paint;
	int fx;
	int fy;
	bool isLayer;
};
Layer layer_identity;
GMatrix matrix_identity;
class EmptyCanvas : public GCanvas {
struct Edge{
	int top_y;
	int bottom_y;
	float curr_x;
	float slope;
	int w;
};

protected:
            virtual void onSaveLayer(const GRect* bounds, const GPaint& paint){ 
				Layer last; //layer before this one that is being created
				///copy bitmap, transform by ctm, calculate new bounds
				last = lastLayer(layerStack);
				GBitmap bm = *new GBitmap();
				//GBitmap bm = last.bounds;
				GRect rect2;
				Layer layer;
				GRect rect;
				if(bounds){
				rect = *bounds;
				//rect = rectFixer(rect);
				
				//rect2 = rectFixer(rect);
				//GRect canvas;
				//canvas.setXYWH(0, 0, last.bounds.width(), last.bounds.height());
				//if(GRoundToInt(rect2.left()) == GRoundToInt(rect2.right()) || GRoundToInt(rect2.top()) == GRoundToInt(rect2.bottom()) || !rect2.intersects(canvas)){
								//return;
				//	}
				//rect2 = rectFixer(rect2);
				//layerRect(layerStack, rect);
			//	Layer last = lastLayer(layers);
				GRect canvas;
				canvas.setXYWH(0, 0, last.bounds.width(), last.bounds.height());
				//rect2 = rectFixer(rect2);
				 if(rect.intersects(canvas)){
				 //setLTRB(T l, T t, T r, T b)
				 rect = rectFixer(rect);
				 //	return;
				 }else{
				 	rect.setLTRB(0, 0, 0, 0);	
				 //	return; 
				 }
				//int arraySize = rect.width()*rect.height();
				//GPixel* pixels[arraySize];
				//GPixel *pointer;
				//pointer = pixels;
				//size_t rb = 0;
				//bm = GBitmap((int) rect.width(), (int) rect.height(), (int) (4*rect.width()), pixels[0], false);
				//bm.reset();
				layer.bounds.alloc(rect.width(), rect.height());
				//if(layer.bounds.width() == 0 && layer.bounds.height() == 0){
				//	layer.bounds = last.bounds;
				//	last.bounds = last.bounds;
				//	layer.fx = last.fx;
				//	layer.fy = last.fy;
				//}else{
				//last.bounds = bm;
				//}
				layer.fx = rect.left();
				layer.fy = rect.top();
				layer.ctm = layerStack.back().ctm;
				layer.ctm.postTranslate(-layer.fx, -layer.fy);
				}else{
					//last.bounds = last.bounds;
					//*bm = last.bounds;
					//layer.bounds = last.bounds;
					layer.bounds.alloc(last.bounds.width(), last.bounds.height());
					layer.ctm = last.ctm;
					//layer.ctm.postTranslate(-last.fx, -last.fy);
					layer.fx = 0;
					layer.fy = 0;
					//rect = GRect::MakeXYWH(last.fx, last.fy, last.bounds.width(), last.bounds.height());
				}
				
				//bm.fWidth =rect2.width();
				//bm.fHeight = rect.height();
				//bm.rb = 0;
				//bm.fPixels = *pixels;
				//bm.isOpaque = false;
				//layer.ctm = layerStack.back().ctm;
				//layer.ctm.postTranslate(-)
				//ctm of new layer = previous then post translate ctm by translation matrix (-L, -T) ctm.postTranslate(-L, -T)
				layer.paint = paint;
				layer.isLayer = true;
				
				layerStack.push_back(layer);
				//onLayerStack.push_back(layer);	
				//onSaveLayer(bounds, paint);
				
			}
public:
    EmptyCanvas(const GBitmap& device) : fDevice(device) {
		
		//last.bounds = this->fDevice;
		matrix_identity.setIdentity();
		layer_identity.bounds = this->fDevice;
		//onLayerStack.push_back(this->fDevice);
		layer_identity.ctm = matrix_identity;
		layer_identity.fx = 0;
		layer_identity.fy = 0;
		layer_identity.isLayer = true;
		layerStack.push_back(layer_identity);
	}

static void blit(int x0, int y, int count, const GPaint& paint, std::vector<Layer>& layerStack){
	//Layer last;
	//
	/*if(x0 < 0){
		x0 = 0;
	}
	*/
	if(x0 + count < x0 || abs(x0) > 10000 || abs(x0 + count) > 10000 || count < 0){
		return;
	}
	Layer last;
	Layer current;
	for(std::vector<Layer>::reverse_iterator i = layerStack.rbegin(); i != layerStack.rend(); ++i){
			current = *i;
			if(current.isLayer){
				last = current;
				//canvas.setXYWH(0, 0, last.bounds.width(), last.bounds.height());
				break;
		//GBitmap(int w, int h, size_t rb, GPixel* pixels, bool isOpaque)
		}
	}
	/*if(count > last.bounds.width() - x0){
		count = last.bounds.width() - x0;
	}
	*/
	//
	//last = lastLayer(layerStack);
	GBlendMode blendMode = paint.getBlendMode();	
	if(paint.getShader()){
			if(!paint.getShader()->setContext(layerStack.back().ctm)){
				return;
			}
		}
	//pipeline: ctm/coverage then shader  	
	if(paint.getShader()){
		//x0 = GRoundToInt(blitLeft.curr_x);
		//x1 = GRoundToInt(blitRight.curr_x);
		//int preCount = x1-x0;
		GPixel storage[count];
		paint.getShader()->shadeRow(x0, y, count, storage);
		for(int j = 0; j < count; j++){
		GPixel *preAddr = last.bounds.getAddr(x0+j, y);
		GPixel prePixel = storage[j];
		/*if(paint.getFilter()){
			paint.getFilter()->filter(&prePixel, &prePixel, 1);
		}
		*/
		GPixel prePixel2 = *preAddr;
		GPixel preBlend = pixelBlender(prePixel, prePixel2, blendMode);
		*preAddr = preBlend;
		}
			
		}else{
		//FILTERTIME
	
		 for(int x = x0; x < x0 + count; x++){
		 	//if(x == clip.right()){
			//	break;
			//}
			
		 	GPixel *addr = last.bounds.getAddr(x, y);
			GPixel pixel = pixelConverter(paint.getColor());
			/*if(paint.getFilter()){
				paint.getFilter()->filter(&pixel, &pixel, 1);
			}
			*/
			if(addr == NULL){
			*addr = pixel;
			}else{
			GPixel pixel2 = *addr;  
			GPixel blend = pixelBlender(pixel, pixel2, blendMode);
					*addr = blend;
		}
		
	}
		
	//place filter
	}
	
	 // end old for loop
	//shadeRow
}


GPoint pQuad(GPoint a, GPoint b, GPoint c, float t){
	float x = (a.x()*(1-t)*(1-t)) + (2*b.x()*t*(1-t))+(c.x()*t*t);
	float y = (a.y()*(1-t)*(1-t)) + (2*b.y()*t*(1-t))+(c.y()*t*t);
	return GPoint::Make(x,y);
}

GPoint pCubic(GPoint a, GPoint b, GPoint c, GPoint d, float t){
	float x = (a.x()*((1-t)*(1-t)*(1-t)))+(3*b.x()*(t*(1-t)*(1-t)))+(3*c.x()*(t*t)*(1-t))+(d.x()*t*t*t);
	float y = (a.y()*((1-t)*(1-t)*(1-t)))+(3*b.y()*(t*(1-t)*(1-t)))+(3*c.y()*(t*t)*(1-t))+(d.y()*t*t*t);
	return GPoint::Make(x,y);
}



virtual void drawPath(const GPath& path, const GPaint& paint) override{
	Layer last;
	last = lastLayer(layerStack);
	GBlendMode blendMode = paint.getBlendMode();
	if(paint.getShader()){
			if(!paint.getShader()->setContext(layerStack.back().ctm)){
				return;
			}
		}
	/*
	struct Edge{
	int top_y;
	int bottom_y;
	float curr_x;
	float slope;
	int w;
};
	*/
	std::vector <Edge> edges;
	Edge e1;
	Edge e2;
	GRect clip;
	clip.setXYWH(0, 0, last.bounds.width(), last.bounds.height());
	GPoint p0;
	GPoint p1;
	
	
	//sort
	GPath p = path;
	p.transform(layerStack.back().ctm);
	GPoint points[4];
	GPoint newPoints[4];
	GPath::Iter iter = GPath::Iter(p);
	GPath::Verb iterCheck = iter.next(points);
	if(iterCheck == GPath::Verb::kDone){
		return;
	}
	GPath::Edger edger = GPath::Edger(p);
	GPath::Verb verbCheck = edger.next(points);
	if(verbCheck == GPath::Verb::kDone){
		return;
	}
 	while(verbCheck != GPath::Verb::kDone){
		if(verbCheck == GPath::Verb::kMove){
			continue;
		}
		if(verbCheck == GPath::Verb::kLine){
//layerStack.back().ctm.mapPoints(newPoints, points, 2);
	clipper(points[0], points[1], edges, clip);
	
	//verbCheck = edger.next(points);
	//continue;
	}
	//verbCheck = edger.next(points);
	if(verbCheck == GPath::Verb::kQuad){
	//clipper(points[0], points[1], edges, clip);
	//clipper(points[1],points[2],edges, clip);
	//verbCheck = edger.next(points);
	//continue;
	
		//#if 0
		GPoint a = points[0];
		GPoint b = points[1];
		GPoint c = points[2];
		GPoint l = GPoint::Make((a.x() - 2*b.x() + c.x())/4.0, (a.y() - 2*b.y() + c.y())/4.0);
		float d = sqrt((l.x()*l.x()) + (l.y() * l.y()));
		int index = ceil(sqrt(4*d));
		float dt = 1.0/index;
		float t = dt;
		/*
		
		
		*/
		
		for(int i = 0; i < index; i++){
			GPoint p0 = pQuad(a, b, c, t);
			GPoint p1 = pQuad(a, b, c, t+dt);
			clipper(p0, p1, edges, clip);
			t+=dt;
			
		}
		//verbCheck = edger.next(points);
		//#endif
	}
	
	if(verbCheck == GPath::Verb::kCubic){
	//clipper(points[0], points[1], edges, clip);
	//clipper(points[1],points[2],edges, clip);
	//clipper(points[2], points[3], edges, clip);
	//verbCheck = edger.next(points);
	//continue;
	GPoint a = points[0];
	GPoint b = points[1];
	GPoint c = points[2];
	GPoint d = points[3];
		GPoint l1 = GPoint::Make(a.x()-2*b.x()+c.x(),a.y()-2*b.y()+c.y());
		GPoint l2 = GPoint::Make(b.x() - 2*c.x() + d.x(), b.y() - 2*c.y() + d.y());
		float d1 = sqrt(l1.x()*l1.x() + l1.y()*l1.y());
		float d2 = sqrt(l2.x()*l2.x() + l2.y()*l2.y());
		GPoint max;
		if(d1 > d2){
			max = l1;
		}else{
			max  = l2;
		}
		float d3 = sqrt(max.x()*max.x() + max.y() * max.y());
		int index = ceil(sqrt((3.0*d3)/(4.0*0.25)));
		float dt = 1.0/index;
		float t = dt;
		for(int i = 0; i < index; i++){
			GPoint p0 = pCubic(a,b,c,d,t);
			GPoint p1 = pCubic(a,b,c,d,t+dt);
			clipper(p0, p1, edges, clip);
			t+=dt;
		}
		
	}
	verbCheck = edger.next(points);
	}
	std::sort (edges.begin(), edges.end(), edgeSort);
	/*for(int i = 0; i < edges.size();  i++){
		if(abs(edges[i].top_y) >= 10000 || abs(edges[0].bottom_y) >=10000 || abs(edges[0].slope >= 10000) || abs(edges[0].curr_x) >= 10000 || abs(edges[0].w) >= 10000){
			for(int j = i; j < edges.size()-1; j++){
					edges[j] = edges[j+1];
				}
				edges.pop_back();
				if(edges.size() <= 1){
					return;
				}
		}
	}
	*/
	
	//WALK AND BLIT
	if(edges.size() < 2){
	
		return;
	}else{
	int x0 = 0;
	int x1 = 0;
	Edge next;
	//Edge edge;
	int w = 0;
	int index = 0;
	for(int y = 0; y < clip.bottom();y++){
	//	if(edges.size() <=1){
		//	return;
		//}
		
		index = 0;
		w = 0;
		//edge = edges[0];
		//x0 = 0;
		//x1 = 0;
		while(index < edges.size() && edges[index].top_y <= y){
			//w+=edge.w;
			if(w == 0){
				x0 = GRoundToInt(edges[index].curr_x);
			}
			//w += edge.w;
			w += edges[index].w;
			if(w == 0){
				x1 = GRoundToInt(edges[index].curr_x);
				if(x0 < 0){
					x0 = 0;
				}
				if(x1 > last.bounds.width()-1){
					x1 = last.bounds.width()-1;
				}
				blit(x0, y, x1-x0, paint, layerStack);
			}
			//next = edges[index+1];
			if(edges[index].bottom_y ==  y+1){
				//erase edge
				for(int i = index; i < edges.size() ; i++){
					edges[i] = edges[i+1];
				}
				edges.pop_back();
		//	if(edges.size() <= 1){
		//			return;
		//		}
		//	edges.erase(edges.begin() + index);
			//	edge = edges[index];
				//index--;
			//	resort_backwards(edges, index);
				//index++;
			}else{
				//edge.curr_x += edge.slope;
				//edges[index] = edge;
				edges[index].curr_x += edges[index].slope;
				index++;
				//resort_backwards(edges, index);
				
			}
				
		//resort_backwards(edges, index);
		
		//edge = edges[index];
		//w+=edge.w;
			
			
		}
		//index = 0;
		//y += 1;
		while(index < edges.size() && edges[index].top_y == y+1){
			//if(edges.size() <= 1){
			//	return;
			//}
			//if(index < edges.size()-1){
			//resort_backwards(edges, index);
				index++;
				//edge = edges[index];
				//
			//}else{
			//	break;
			//}
			
			//next = edges[index+1];
			
			
		}
		resort_backwards(edges, index);
		
		}
		 
}        
}



void save() override{
	Layer layer;
	layer.ctm = layerStack.back().ctm;
	layer.isLayer = false;
	layerStack.push_back(layer);
}
void restore() override{ 

Layer back = layerStack.back();
	if(layerStack.size() >= 1){
	layerStack.pop_back();
	}
	if(back.isLayer){
	Layer last = lastLayer(layerStack);
	last.ctm.postTranslate(back.fx, back.fy);
	if(back.paint.getShader()){
		if(!back.paint.getShader()->setContext(last.ctm)){
			return;
		}
	}
	
		GRect deviceRect;
	//deviceRect.setXYWH(back.fx, back.fy, back.bounds.width(), back.bounds.height());
		
	
		//GBitmap back.bounds = back.bounds;
		if(back.paint.getShader()){
			for (int i = 0; i < back.bounds.width(); i++){
			for(int j = 0; j < back.bounds.height(); j++){
				GPixel *src = back.bounds.getAddr(i, j);
				GPixel *dst = last.bounds.getAddr(i + back.fx, j + back.fy);
				/**
     *  Given a row of pixels in device space [x, y] ... [x + count - 1, y], return the
     *  corresponding src pixels in row[0...count - 1]. The caller must ensure that row[]
     *  can hold at least [count] entries.
     */
				//GPixel blend = pixelBlender(*src, *dst, back.paint.getBlendMode());
				GPixel* storage;
				back.paint.getShader()->shadeRow(i+ back.fy, j +back.fy, 1, storage);
				GPixel blend;
				blend = *storage;
				GPixel filtered;
				//if(back.paint.getFilter()){
				//check filter
			//	back.paint.getFilter()->filter(&blend, &blend, 1);
				//*dst = blend;
			//	}
				GPixel blend2 = pixelBlender(blend, *dst, back.paint.getBlendMode());
				*dst = blend2;
			}
		}
		}else{
		//have to draw rect manually since there are two bitmaps to consider when drawing
		for (int i = 0; i < back.bounds.width(); i++){
			for(int j = 0; j < back.bounds.height(); j++){
				GPixel *src = back.bounds.getAddr(i, j);
				GPixel *dst = last.bounds.getAddr(i + back.fx, j + back.fy);
				//if(GPixel_GetA(*dst) == 0){
				//GPixel GPixel_PackARGB(unsigned a, unsigned r, unsigned g, unsigned b)
				//	*dst = GPixel_PackARGB(0, 0, 0, 0);
				//}
				
				GPixel filtered = *src;
				//if(back.paint.getFilter()){
				//back.paint.getFilter()->filter(&filtered, &filtered, 1);
				//*dst = blend;
				//}
				GPixel blend = pixelBlender(filtered, *dst, back.paint.getBlendMode());
				*dst = blend;
			}
		}
		}
		//drawRect(deviceRect, back.paint);
			//layerStack = restoreStack;
	//	last.bounds = last.bounds;
		return;
		//onLayerStack.pop_back();
		//onSaveLayer(onLayerStack.back().bounds, onLayerStack.back().paint);
	
	}else{
		return;
	}
	
}
void concat(const GMatrix& matrix) override{
/**
     *  Modifies the CTM by preconcatenating the specified matrix with the CTM. The canvas
     *  is constructed with an identity CTM.
     *
     *  CTM' = CTM * matrix
     */
	GMatrix back = layerStack.back().ctm;
	layerStack.back().ctm.setConcat(back, matrix);
}


    void drawPaint(const GPaint& paint) override{ 
	if(paint.getShader()){
		if(!paint.getShader()->setContext(layerStack.back().ctm)){
			return;
		}
	}
	Layer last;
	last = lastLayer(layerStack);
	GBlendMode blendMode = paint.getBlendMode();
	//consider filter then blend
	GPixel* pixels[last.bounds.width()];
	GPixel output;
	int preCount = last.bounds.width();
	if(paint.getShader()){
		GPixel storage[preCount];
			/**
     *  Given a row of pixels in device space [x, y] ... [x + count - 1, y], return the
     *  corresponding src pixels in row[0...count - 1]. The caller must ensure that row[]
     *  can hold at least [count] entries.
     */
	// 	for(int i = 0; i < last.bounds.height(); i++){
		
	//	}
		
		//if(paint.getFilter()){
		//for(int i = 0; i < last.bounds.height(); i++){
		//for(int j = 0; j < last.bounds.width(); j++){
		//	pixels[j] = last.bounds.getAddr(j, i);
		//	}
		//}
		//}
		for(int i = 0; i < last.bounds.height(); i++){
		paint.getShader()->shadeRow(0, i, preCount, storage);
		for(int j = 0; j < last.bounds.width(); j++){
		GPixel *preAddr = last.bounds.getAddr(j, i);
		GPixel prePixel = storage[j];
		GPixel prePixel2 = *preAddr;
		/*if(paint.getFilter()){
			paint.getFilter()->filter(&prePixel, &prePixel, 1);
		}
		*/
		GPixel preBlend = pixelBlender(prePixel, prePixel2, blendMode);
		*preAddr = preBlend;
		}
		}
		
			
			 //for(int x = 0; x < last.bounds.width(); x++){
	  		 	
		//}
	
	}
	else{	
      for(int x = 0; x < last.bounds.width(); x++){
	  	for(int y = 0; y < last.bounds.height(); y++){
		GPixel *addr = last.bounds.getAddr(x, y);
			GPixel pixel = pixelConverter(paint.getColor());
			/*if(paint.getFilter()){
				paint.getFilter()->filter(&pixel, &pixel, 1);
			}
			*/
			if(addr == NULL){
			*addr = pixel;
			}else{
			GPixel pixel2 = *addr;
			*addr = pixelBlender(pixel, pixel2, blendMode);
			}
		}
	 }
	 
	 }
	 
		
    }
    //SetXYWH- left -x, top -  y, right -x+w, bottom - y + h
	//SetLTRB - left l, top t, right r, bottom b
	
    void drawRect(const GRect& rect, const GPaint& paint) override {
        // your code here
		Layer last;
		last = lastLayer(layerStack);
		GBlendMode blendMode;
		GRect canvas;
		canvas.setXYWH(0, 0, last.bounds.width(), last.bounds.height());
		GPoint fourPoints[4];
		GPoint p1;
		GPoint p2;
		GPoint p3;
		GPoint p4;
		if(!rect.intersects(canvas)){
			return;
		}
		//rect.left and rect.right +x offset, rec.top andrect.bottom + yoffset
		GRect rect2;
		rect2.setLTRB(GRoundToInt(rect.left()), GRoundToInt(rect.top()), GRoundToInt(rect.right()), GRoundToInt(rect.bottom()));
		if(GRoundToInt(rect2.left()) == GRoundToInt(rect2.right()) || GRoundToInt(rect2.top()) == GRoundToInt(rect2.bottom()) || !rect2.intersects(canvas)){
			return;
		}
		rect2 = rectFixer(rect);
		p1.set(rect2.left(), rect2.top());
		p2.set(rect2.right(), rect2.top());
		p3.set(rect2.right(), rect2.bottom());
		p4.set(rect2.left(), rect2.bottom());
		
		fourPoints[0] = p1;
		fourPoints[1] = p2;
		fourPoints[2] = p3;
		fourPoints[3] = p4;
		drawConvexPolygon(fourPoints, 4, paint);

}
//if a.top != b.top  ret a.top < b.top .... curr_x, slope
//static inline bool pathSor(const Edge& e1, Edge& e2){
	
//}

//begin methods
class TriColorShader : public GShader {
	public:
		TriColorShader(GPoint p0, GPoint p1, GPoint p2, const GColor colors[]) : po0(p0), po1(p1), po2(p2){
		sColors = new GColor[3]; 
		memcpy(sColors, colors, 3*sizeof(GColor));
		float ux = po1.x() - po0.x();
		float uy = po1.y() - po0.y();
		float vx = po2.x() - po0.x();
		float vy = po2.y() - po0.y();
		lM.set6(ux, vx, po0.x(), uy, vy, po0.y());
		}
	private:
		GPoint po0;
		GPoint po1;
		GPoint po2;
		GColor* sColors;
		GMatrix lM;
		GMatrix current;
	
	bool isOpaque() override{
		for(int i = 0; i < 3; i++){
		GPixel curr = pixelConverter(sColors[i]);
		if(GPixel_GetA(curr) == 1){
			continue;
		}
		return false;
		}
		return true;
	}
	/*
	GMatrix lInv;
		localMatrix.invert(&lInv);
		// lInv.invert(&localMatrix);
		//temp.setConcat(ctm, localMatrix);
		if(!ctm.invert(&current)){
			return false;
		}
		current.postConcat(lInv);
	*/
	bool setContext(const GMatrix& ctm) override{
		GMatrix lInv;
		GMatrix c;
		c.setConcat(ctm, lM);
		if(!c.invert(&lInv)){
		return false;
		}
		current = lInv;
		return true;
	}
	 void shadeRow(int x, int y, int count, GPixel row[]) override{
	 	 for(int i = 0; i < count; i++){
			 	GPoint localPoint = current.mapXY(x + 0.5 + i, y + 0.5);
				GColor gradient;
				GColor c0 = sColors[0];
				GColor c1 = sColors[2];
				GColor c2 = sColors[1];
				float a = (localPoint.x()*c2.fA) + (localPoint.y()*c1.fA) + ((1-localPoint.x() - localPoint.y())*c0.fA);
				float r =  (localPoint.x()*c2.fR) + (localPoint.y()*c1.fR) + ((1-localPoint.x() - localPoint.y())*c0.fR);
				float g =  (localPoint.x()*c2.fG) + (localPoint.y()*c1.fG) + ((1-localPoint.x() - localPoint.y())*c0.fG);
				float b =  (localPoint.x()*c2.fB) + (localPoint.y()*c1.fB) + ((1-localPoint.x() - localPoint.y())*c0.fB);
				gradient.fA = a;
				gradient.fR = r;
				gradient.fG = g;
				gradient.fB = b;
				row[i] = pixelConverter(gradient);
	 }		
  }
	  static inline GPixel pixelConverter(const GColor color){
			//GColor color = paint.getColor();
			GColor pinned = color.pinToUnit();
			float a = pinned.fA* 255;
			int a2 = GRoundToInt(a);
			float r = pinned.fR*pinned.fA;
			r = r*255;
			int r2 = GRoundToInt(r);
			float g = pinned.fG*pinned.fA;
			g = g*255;
			int g2 = GRoundToInt(g);
			float b = pinned.fB*pinned.fA;
			b = b*255;
			int b2 = GRoundToInt(b);
			return GPixel_PackARGB(a2, r2, g2, b2);
		} 	
  
};
GMatrix compute_basis(GPoint pt0, GPoint pt1, GPoint pt2){
	GMatrix m;
	float ux = pt1.x() - pt0.x();
	float uy = pt1.y() - pt0.y();
	float vx = pt2.x() - pt0.x();
	float vy = pt2.y() - pt0.y();
	m.set6(ux, vx, pt0.x(), uy, vy, pt0.y());
	return m;
}
void drawTriangleWithTex(const GPoint pts[3], const GPoint tex[3], const GPaint& paint){
	GMatrix p;
	GMatrix t;
	GMatrix invT;
	p = compute_basis(pts[0], pts[1], pts[2]);
	t = compute_basis(tex[0], tex[1], tex[2]);
	t.invert(&invT);
	//"something"
	GMatrix mp;
	mp.setConcat(p, invT);
	ProxyShader proxy(paint.getShader(), mp);
	drawConvexPolygon(pts, 3, GPaint(&proxy));
}
class ProxyShader : public GShader {
	public:
		ProxyShader(GShader* shader, const GMatrix& extraTransform) : fRealShader(shader), fExtraTransform(extraTransform) {}
	private:
		GShader* fRealShader;
		GMatrix fExtraTransform;
	bool isOpaque() override {return fRealShader->isOpaque();}
	bool setContext(const GMatrix& ctm) override {
		GMatrix ctm_extra;
		ctm_extra.setConcat(ctm, fExtraTransform);
		return fRealShader->setContext(ctm_extra);
	}
	void shadeRow(int x, int y, int count, GPixel row[]) override {
		fRealShader->shadeRow(x, y, count, row);
	}	
};
class TriTexShader : public GShader {
	public:
		TriTexShader(const GBitmap& device, GPoint p0, GPoint p1, GPoint p2, const GPoint points[]) : sDevice(device), po0(p0), po1(po1), po2(po2){
		sPoints = new GPoint[3]; 
		memcpy(sPoints, points, 3*sizeof(GPoint));
		}
	private:
		GPoint po0;
		GPoint po1;
		GPoint po2;
		GPoint* sPoints;
		const GBitmap& sDevice;
		GMatrix current;
	
	bool isOpaque() override{
		return true;
	}
	bool setContext(const GMatrix& ctm) override{
		GMatrix lM;
		GPoint t0 = sPoints[0];
		GPoint t1 = sPoints[1];
		GPoint t2 = sPoints[2];
		float t1x = t1.x() - t0.x();
		float t1y = t1.y() - t0.y();
		float t2x = t2.x() - t0.x();
		float t2y = t2.y() - t0.y();
		lM.set6(t1x, t2x, t0.x(), t1y, t2y, t0.y());
		GMatrix p;
		//float total_x = po0.x() + po1.x() + po2.x();
		//float total_y = po0.y() +  po1.y() + po2.y();
		float ux = po1.x() - po0.x();
		float uy = po1.y() - po0.y();
		float vx = po2.x() - po0.x();
		float vy = po2.y() - po0.y();
		p.set6(ux, vx, po0.x(), uy, vy, po0.y());
		GMatrix lInv;
		lM.invert(&lInv);
		GMatrix pc;
		pc.setConcat(p, lInv);
		// lInv.invert(&localMatrix);
		//temp.setConcat(ctm, localMatrix);
		current.setConcat(ctm, pc);
		return true;	
	}
	 void shadeRow(int x, int y, int count, GPixel row[]) override{
	 	 for(int i = 0; i < count; i++){
			 	GPoint localPoint = current.mapXY(x + 0.5 +  i, y+ 0.5);
				GPixel* pixel = this->sDevice.getAddr(x, y);
				row[i] = *pixel;
	 }		
  }
};
class TriComposeShader : public GShader {
	public:
		TriComposeShader(GShader* sh1, GShader* sh2) : shader1(sh1), shader2(sh2) {}
	private:
	//void shadeRow(int x, int y, int count, GPixel row[]) override{	
		GShader* shader1;
		GShader* shader2;
		GMatrix current;
	
	bool isOpaque() override{
		return true;
	}
		
	bool setContext(const GMatrix& ctm) override{
		if(shader1->setContext(ctm) && shader2->setContext(ctm)){
			return true;
		}else{
			return false;
		}
	}
	
	void shadeRow(int x, int y, int count, GPixel row[]) override{
		GPixel* storage1;
		GPixel* storage2;
		storage1 = new GPixel[count];
		storage2 = new GPixel[count];
		memcpy(storage1, row, count*sizeof(GPixel));
		memcpy(storage2, row, count*sizeof(GPixel));
		//GPoint* sPoints;
		//sColors = new GColor[3]; 
		//memcpy(sColors, colors, 3*sizeof(GColor));
		shader1->shadeRow(x, y, count, storage1);
		shader2->shadeRow(x, y, count, storage2);
		for(int i = 0; i < count; i++){
			float a = GRoundToInt((GPixel_GetA(storage1[i]) * GPixel_GetA(storage2[i]))/255.0f);
			float r = GRoundToInt((GPixel_GetR(storage1[i]) * GPixel_GetR(storage2[i]))/255.0f);
			float g = GRoundToInt((GPixel_GetG(storage1[i]) * GPixel_GetG(storage2[i]))/255.0f);
			float b = GRoundToInt((GPixel_GetB(storage1[i]) * GPixel_GetB(storage2[i]))/255.0f);
			row[i] = GPixel_PackARGB(a, r, g, b);
		}
		
	}
};
	 static GShader* GCreateTriColorShader(GPoint p0, GPoint p1, GPoint p2, const GColor colors[3]){
		return new TriColorShader(p0, p1, p2, colors);
		//device w/localInv
	}
	 static GShader* GCreateTriTexShader(const GBitmap& device, GPoint p0, GPoint p1, GPoint p2, const GPoint points[3]){
			if(!device.pixels()){
			return nullptr;
		}
		return new TriTexShader(device, p0, p1, p2, points);
		//device w/localInv
	}
	 static GShader* GCreateTriComposeShader(GShader* sh1, GShader* sh2){
		return new TriComposeShader(sh1, sh2);
		//device w/localInv
	}
//end methods 11/26
     /**
         *  Draw a mesh of triangles, with optional colors and/or texture-coordinates at each vertex.
         *
         *  The triangles are specified by successive triples of indices.
         *      int n = 0;
         *      for (i = 0; i < count; ++i) {
         *          point0 = vertx[indices[n+0]]
         *          point1 = verts[indices[n+1]]
         *          point2 = verts[indices[n+2]]
         *          ...
         *          n += 3
         *      }
         *
         *  If colors is not null, then each vertex has an associated color, to be interpolated
         *  across the triangle. The colors are referenced in the same way as the verts.
         *          color0 = colors[indices[n+0]]
         *          color1 = colors[indices[n+1]]
         *          color2 = colors[indices[n+2]]
         *
         *  If texs is not null, then each vertex has an associated texture coordinate, to be used
         *  to specify a coordinate in the paint's shader's space. If there is no shader on the
         *  paint, then texs[] should be ignored. It is referenced in the same way as verts and colors.
         *          texs0 = texs[indices[n+0]]
         *          texs1 = texs[indices[n+1]]
         *          texs2 = texs[indices[n+2]]
         *
         *  If both colors and texs[] are specified, then at each pixel their values are multiplied
         *  together, component by component.
         */
//  virtual void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[],
//                          int count, const int indices[], const GPaint& paint) = 0;
 //void drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) override 
virtual void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count, const int indices[], const GPaint& paint) override {
			    int n = 0;
			    //local matrix with triangle points
				//u v p0
				Layer last;
				last = lastLayer(layerStack);
			    GBlendMode blendMode = paint.getBlendMode();
			    GColor color0;
				GColor color1;
				GColor color2;
				GPoint texs0;
				GPoint texs1;
				GPoint texs2;
	 			bool colors_bool = false;
				bool texs_bool = false;
				if(colors != NULL){
					colors_bool = true;
				}
				
				if(texs != NULL){
					texs_bool = true;
				}
				
               for (int i = 0; i < count; ++i) {
	                   GPoint point0 = verts[indices[n+0]];
	                   GPoint  point1 = verts[indices[n+1]];
	                   GPoint point2 = verts[indices[n+2]];
				  if(colors_bool){
					   GColor color0 = colors[indices[n+0]];
	                   GColor color1 = colors[indices[n+1]];
	                   GColor color2 = colors[indices[n+2]];
				   }
				  if(texs_bool){
					   GPoint texs0 = texs[indices[n+0]];
	                   GPoint texs1 = texs[indices[n+1]];
	                   GPoint texs2 = texs[indices[n+2]];
				  }
				   GPoint curr_points[3];
				   curr_points[0] = verts[indices[n+0]];
				   curr_points[1] = verts[indices[n+1]];
				   curr_points[2] = verts[indices[n+2]];
				   GPoint curr_texs[3];
				   if(texs_bool){
				   curr_texs[0] = texs[indices[n+0]];
				   curr_texs[1] = texs[indices[n+1]];
				   curr_texs[2] = texs[indices[n+2]];
				   }
				   GColor curr_colors[3];
				   if(colors_bool){
				   curr_colors[0] = colors[indices[n+0]];
				   curr_colors[1] = colors[indices[n+1]];
				   curr_colors[2] = colors[indices[n+2]];
				   }
				 //drawConvexPolygon(curr, 3, paint);
				 //void drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) override 
				// auto tex;
				 //auto color;
				 //auto c1;
				 //auto c2;
				 GPaint shadedPaint = paint;
				 if(colors_bool && texs_bool){
				 //points array holds textures
				 //11 26
				 //GCreateTriComposeShader
				 	GMatrix p;
					GMatrix t;
					GMatrix invT;
					p = compute_basis(verts[indices[n+0]], verts[indices[n+1]], verts[indices[n+2]]);
					t = compute_basis(texs[indices[n+0]],  texs[indices[n+1]], texs[indices[n+2]]);
					t.invert(&invT);
					//"something"
					GMatrix mp;
					mp.setConcat(p, invT);
					ProxyShader proxy(paint.getShader(), mp);
					//auto tex = GCreateTriTexShader(last.bounds, point0, point1, point2, curr_texs);
					TriColorShader colsh(point0, point1, point2, curr_colors);
					TriComposeShader compose(&proxy, &colsh);
				 	shadedPaint.setShader(&compose);
					//if(!shadedPaint.getShader()->setContext(last.ctm)){
					//			continue;
					//		}
				 	drawConvexPolygon(curr_points, 3, shadedPaint);
				}
				else if(colors_bool && !texs_bool){
					TriColorShader colsh(point0, point1, point2, curr_colors);
					/*if(paint.getShader()){
				 	auto c1 = GCreateTriComposeShader(color, paint.getShader());	
					shadedPaint.setShader(c1);
					//if(!shadedPaint.getShader()->setContext(last.ctm)){
					//			continue;
					//		}
					drawConvexPolygon(curr_points, 3, shadedPaint);
				 }else{
				 */
				 	//shadedPaint.setColor(color0);
					shadedPaint.setShader(&colsh);
					//if(!shadedPaint.getShader()->setContext(last.ctm)){
					//			continue;
					//		}
					
					drawConvexPolygon(curr_points, 3, shadedPaint);
					
				}else if(!colors_bool && texs_bool){
					drawTriangleWithTex(curr_points, curr_texs, paint);
				//}else{
				}	
				//}			  
                   n += 3;
               }
	
}

// virtual void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
//                          int level, const GPaint&) = 0;
 /**
         *  Draw the quad, with optional color and/or texture coordinate at each corner. Tesselate
         *  the quad based on "level":
         *      level == 0 --> 1 quad  -->  2 triangles
         *      level == 1 --> 4 quads -->  8 triangles
         *      level == 2 --> 9 quads --> 18 triangles
         *      ...
         *  The 4 corners of the quad are specified in this order:
         *      top-left --> top-right --> bottom-right --> bottom-left
         *  Each quad is triangulated on the diagonal top-right --> bottom-left
         *      0---1
         *      |  /|
         *      | / |
         *      |/  |
         *      3---2
         *
         *  colors and/or texs can be null. The resulting triangles should be passed to drawMesh(...).
         */
virtual void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4], int level, const GPaint& paint) override{
				bool colors_bool = false;
				bool texs_bool = false;
				if(colors != NULL){
					colors_bool = true;
				}
				
				if(texs != NULL){
					texs_bool = true;
				}
	GPoint a = verts[0];
	GPoint b = verts[1];
	GPoint c = verts[2];
	GPoint d = verts[3];
	GColor e;
	GColor f;
	GColor g;
	GColor h;
	GPoint i_p;
	GPoint j_p;
	GPoint k_p;
	GPoint l_p;
	if(colors_bool){
		e = colors[0];
		f = colors[1];
		g = colors[2];
		h = colors[3];
	}
	if(texs_bool){
		i_p = texs[0];
		j_p = texs[1];
		k_p = texs[2];
		l_p = texs[3];
	}
	/*GMatrix m0;
	float width = b.x();
	float height = d.y();
	m0.set6(1/width, 0, 0, 0, 1/height, 0);
	GMatrix m1;
	m1.set6(b.x()-a.x(), c.x()-a.x(), a.x(), b.y()-a.y(), c.y()-a.y(), a.y());
	GMatrix dm;
	dm.setConcat(m1, m0);
	GMatrix dminv;
	dm.invert(&dminv);
	*/
	GPoint totalVerts[10000];
	GColor totalColors[10000];
	GPoint totalTexs[10000];
	float u = 0;
	float v = 0;
	float increment = 1.0f/(level+1);
	GPoint curr_points[4];
	GColor curr_colors[4];
	GPoint curr_texs[4];
	int k = 0;
	for(int i = 0; i <= level+1; i++){
		v = 0;
		for(int j = 0; j < level+1; j++){
			//totalVerts[k] = GPoint::Make((1-v)*(((1-u)*a.x())+(u*b.x())) + v*(((1-u)*d.x()) + (u*c.x())),(1-v)*(((1-u)*a.y())+(u*b.y())) + v*(((1-u)*d.y()) + (u*c.y())));
			curr_points[0] = GPoint::Make((1-v)*(((1-u)*a.x())+(u*b.x())) + v*(((1-u)*d.x()) + (u*c.x())),(1-v)*(((1-u)*a.y())+(u*b.y())) + v*(((1-u)*d.y()) + (u*c.y())));
			curr_points[1] = GPoint::Make((1-(v+increment))*(((1-u)*a.x())+(u*b.x())) + (v+increment)*(((1-u)*d.x()) + (u*c.x())),(1-(v+increment))*(((1-u)*a.y())+(u*b.y())) + (v+increment)*(((1-u)*d.y()) + (u*c.y())));
			curr_points[2] = GPoint::Make((1-v)*(((1-(u+increment))*a.x())+((u+increment)*b.x())) + v*(((1-(u+increment))*d.x()) + ((u+increment)*c.x())),(1-v)*(((1-(u+increment))*a.y())+((u+increment)*b.y())) + v*(((1-(u+increment))*d.y()) + ((u+increment)*c.y())));
			curr_points[3] = GPoint::Make((1-(v+increment))*(((1-(u+increment))*a.x())+((u+increment)*b.x())) + (v+increment)*(((1-(u+increment))*d.x()) + ((u+increment)*c.x())),(1-(v+increment))*(((1-(u+increment))*a.y())+((u+increment)*b.y())) + (v+increment)*(((1-(u+increment))*d.y()) + ((u+increment)*c.y())));
			if(colors_bool){
				GColor c1;
				GColor c2;
				GColor c3;
				GColor c4;
				c1.fA = (1-v)*(((1-u)*e.fA)+(u*f.fA)) + v*(((1-u)*h.fA) + (u*g.fA));
				c1.fR = (1-v)*(((1-u)*e.fR)+(u*f.fR)) + v*(((1-u)*h.fR) + (u*g.fR));
				c1.fG = (1-v)*(((1-u)*e.fG)+(u*f.fG)) + v*(((1-u)*h.fG) + (u*g.fG));
				c1.fB = (1-v)*(((1-u)*e.fB)+(u*f.fB)) + v*(((1-u)*h.fB) + (u*g.fB));
				//totalColors[k] = c1;
				curr_colors[0] = c1;
				c2.fA = (1-(v+increment))*(((1-u)*e.fA)+(u*f.fA)) + (v+increment)*(((1-u)*h.fA) + (u*g.fA));
				c2.fR = (1-(v+increment))*(((1-u)*e.fR)+(u*f.fR)) + (v+increment)*(((1-u)*h.fR) + (u*g.fR));
				c2.fG = (1-(v+increment))*(((1-u)*e.fG)+(u*f.fG)) + (v+increment)*(((1-u)*h.fG) + (u*g.fG));
				c2.fB = (1-(v+increment))*(((1-u)*e.fB)+(u*f.fB)) + (v+increment)*(((1-u)*h.fB) + (u*g.fB));
				curr_colors[1] = c2;
				c3.fA = (1-v)*(((1-(u+increment))*e.fA)+((u+increment)*f.fA)) + v*(((1-(u+increment))*h.fA) + ((u+increment)*g.fA));
				c3.fR = (1-v)*(((1-(u+increment))*e.fR)+((u+increment)*f.fR)) + v*(((1-(u+increment))*h.fR) + ((u+increment)*g.fR));
				c3.fG = (1-v)*(((1-(u+increment))*e.fG)+((u+increment)*f.fG)) + v*(((1-(u+increment))*h.fG) + ((u+increment)*g.fG));
				c3.fB = (1-v)*(((1-(u+increment))*e.fB)+((u+increment)*f.fB)) + v*(((1-(u+increment))*h.fB) + ((u+increment)*g.fB));
				curr_colors[2] = c3;
				c4.fA = (1-(v+increment))*(((1-(u+increment))*e.fA)+((u+increment)*f.fA)) + (v+increment)*(((1-(u+increment))*h.fA) + ((u+increment)*g.fA));
				c4.fR = (1-(v+increment))*(((1-(u+increment))*e.fR)+((u+increment)*f.fR)) + (v+increment)*(((1-(u+increment))*h.fR) + ((u+increment)*g.fR));
				c4.fG = (1-(v+increment))*(((1-(u+increment))*e.fG)+((u+increment)*f.fG)) + (v+increment)*(((1-(u+increment))*h.fG) + ((u+increment)*g.fG));
				c4.fB = (1-(v+increment))*(((1-(u+increment))*e.fB)+((u+increment)*f.fB)) + (v+increment)*(((1-(u+increment))*h.fB) + ((u+increment)*g.fB));
				curr_colors[3] = c4;
			}
			if(texs_bool){
				//totalTexs[k] = GPoint::Make((1-v)*(((1-u)*i_p.x())+(u*j_p.x())) + v*(((1-u)*l_p.x()) + (u*k_p.x())),(1-v)*(((1-u)*i_p.y())+(u*j_p.y())) + v*(((1-u)*l_p.y()) + (u*k_p.y())));
				curr_texs[0] = GPoint::Make((1-v)*(((1-u)*i_p.x())+(u*j_p.x())) + v*(((1-u)*l_p.x()) + (u*k_p.x())),(1-v)*(((1-u)*i_p.y())+(u*j_p.y())) + v*(((1-u)*l_p.y()) + (u*k_p.y())));
				curr_texs[1] = GPoint::Make((1-(v+increment))*(((1-u)*i_p.x())+(u*j_p.x())) + (v+increment)*(((1-u)*l_p.x()) + (u*k_p.x())),(1-(v+increment))*(((1-u)*i_p.y())+(u*j_p.y())) + (v+increment)*(((1-u)*l_p.y()) + (u*k_p.y())));
				curr_texs[2] =  GPoint::Make((1-v)*(((1-(u+increment))*i_p.x())+((u+increment)*j_p.x())) + v*(((1-(u+increment))*l_p.x()) + ((u+increment)*k_p.x())),(1-v)*(((1-(u+increment))*i_p.y())+((u+increment)*j_p.y())) + v*(((1-(u+increment))*l_p.y()) + ((u+increment)*k_p.y())));
				curr_texs[3] =  GPoint::Make((1-(v+increment))*(((1-(u+increment))*i_p.x())+((u+increment)*j_p.x())) + (v+increment)*(((1-(u+increment))*l_p.x()) + ((u+increment)*k_p.x())),(1-(v+increment))*(((1-(u+increment))*i_p.y())+((u+increment)*j_p.y())) + (v+increment)*(((1-(u+increment))*l_p.y()) + ((u+increment)*k_p.y())));
			}
			int indices[6];
			indices[0] = 0;
			indices[1] = 1;
			indices[2] = 2;
			indices[3] = 2;
			indices[4] = 1;
			indices[5] = 3;
		if(colors_bool && texs_bool){
			drawMesh(curr_points, curr_colors, curr_texs, 2, indices, paint);
		}else if(colors_bool && !texs_bool){
			drawMesh(curr_points, curr_colors, nullptr, 2, indices, paint);
		}else if(!colors_bool && texs_bool){
			drawMesh(curr_points, nullptr, curr_texs, 2, indices, paint);
		}else{
			return;
		}
			
			v += increment;
			k++;
		}
		
		u += increment;
	}
	/*int numberofTriangles = 2*(level+1)*(level+1);
	int indices[numberofTriangles*3];
	for(int l = 0; l < numberofTriangles*3; l++){
		indices[l] = l;
	}
	//dminv.mapPoints(totalVerts, numberofTriangles*3);
	if(colors_bool && texs_bool){
		drawMesh(totalVerts, totalColors, totalTexs, numberofTriangles, indices, paint);	
	}else if(colors_bool && !texs_bool){
		drawMesh(totalVerts, totalColors, nullptr, numberofTriangles, indices, paint);	
	}else if(!colors_bool && texs_bool){
		drawMesh(totalVerts, nullptr, totalTexs, numberofTriangles, indices, paint);	
	}else{
		return;
	}
	*/
		//drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count, const int indices[], const GPaint& paint)
		
	//GPoint newVerts[4];
	//dm.mapPoints(verts, newVerts, 4);
	/*for(int i = 0; i < level; i++){
		for(int j = 0; j < 4^i){
		GPoint totalVerts[10000];
		//number of quads = level*4^(i)
		//u/2 for 1st level u/4 for 2nd u/8 ...
		totalVerts[0] = a;
		totalVerts[]	
		}
	}
	*/
	return;
}
void drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) override {
Layer last;
last = lastLayer(layerStack);

GBlendMode blendMode = paint.getBlendMode();
if(paint.getShader()){
		if(!paint.getShader()->setContext(layerStack.back().ctm)){
			return;
		}
	}
	
	std::vector <Edge> edges;
	Edge e1;
	Edge e2;
	//GBlendMode blendMode = paint.getBlendMode()
	GRect clip;
	clip.setXYWH(0, 0, last.bounds.width(), last.bounds.height());
	GPoint p0;
	GPoint p1;
	//void mapPoints(GPoint& dst[], const GPoint src[], int count) const override{
	GPoint newPoints[count];

	layerStack.back().ctm.mapPoints(newPoints, points, count);
	
	for(int i = 0; i < count; i++){
		if(i == count-1){
		p0 = newPoints[i];
		p1 = newPoints[0]; 
		clipper(p0, p1, edges, clip);	
		}else{
			p0 = newPoints[i];
			p1 = newPoints[i+1];
			clipper(p0, p1, edges, clip);	
		}
	}
	
	//sort
	std::sort (edges.begin(), edges.end(), edgeSort);
	
	//WALK AND BLIT
	if(edges.size() < 2){
		//using namespace std;
		//std::printf("edges less than 2");
		return;
	}else{
	Edge blitLeft = edges[0];
	Edge blitRight = edges[1];
	//edge.top_y and edge.bot_y += yoffset
	//edge.curr_x += xoffset
	//edge.slope(edited in makeedge) x1+offset -x0 + offset/y1+offset - y0 +offset
	//edge.bot_y = bot_y*yscale
	//edge.slope = x_bottom_y*xscale-x0/y1*yscale
	int y = blitLeft.top_y;
	int i = 1;
	int x0;
	int x1;
	while(y < clip.bottom()){
		if(GRoundToInt(blitLeft.curr_x) < clip.left()){
			blitLeft.curr_x = 0;
		}
		
	//pipeline: ctm/coverage then shader  	
	if(paint.getShader()){
		x0 = GRoundToInt(blitLeft.curr_x);
		x1 = GRoundToInt(blitRight.curr_x);
		int preCount = x1-x0;
		GPixel storage[preCount];
		paint.getShader()->shadeRow(x0, y, preCount, storage);
		for(int j = 0; j < preCount; j++){
		GPixel *preAddr = last.bounds.getAddr(x0+j, y);
		GPixel prePixel = storage[j];
		/*if(paint.getFilter()){
			paint.getFilter()->filter(&prePixel, &prePixel, 1);
		}
		*/
		GPixel prePixel2 = *preAddr;
		GPixel preBlend = pixelBlender(prePixel, prePixel2, blendMode);
		*preAddr = preBlend;
		}
			
		}else{
		//FILTERTIME
	
		 for(int x = GRoundToInt(blitLeft.curr_x); x < GRoundToInt(blitRight.curr_x); x++){
		 	if(x == clip.right()){
				break;
			}
			
		 	GPixel *addr = last.bounds.getAddr(x, y);
			GPixel pixel = pixelConverter(paint.getColor());
			/*if(paint.getFilter()){
				paint.getFilter()->filter(&pixel, &pixel, 1);
			}
			*/
			if(addr == NULL){
			*addr = pixel;
			}else{
			GPixel pixel2 = *addr;  
			GPixel blend = pixelBlender(pixel, pixel2, blendMode);
					*addr = blend;
		}
		
	}
		
	//place filter
	}
	
	 // end old for loop
	//shadeRow
	
	y = y+1;
	blitLeft.curr_x += blitLeft.slope;
	blitRight.curr_x += blitRight.slope;
	if(y == blitLeft.bottom_y && y == blitRight.bottom_y){
if(i == edges.size()-1){
	return;
}
	blitLeft = edges[i+1];
	blitRight = edges[i+2];
	i = i+2;
	}
	
	else if(y == blitLeft.bottom_y){
		if(i == edges.size()){
			return;
		}
		blitLeft = edges[i+1];
		i = i + 1;
	}
	else if(y == blitRight.bottom_y){
		if(i == edges.size()){
			return;
		}
		blitRight = edges[i+1];
		i = i + 1;
	}
	
}
//outside while loop	
	
	
}
//if left bot <= y left .-< bottom

}	
static inline void resort_backwards(std::vector<Edge>& edges, int& index){
	
	if(index == 0){
		//edges[0].w *= -1;
		return;
	}
	int curr = index;
	//edges[index].w *= -1;
	/*
	int curr = index;
	Edge e = edges[index];
	//std::vector<Edge>:: iterator prv = ptr -1;
	//std::vector<Edge>:: iterator prv = ptr -1;
	Edge prev = edges[index-1]; 
	bool sorted = false;
	while(!sorted){
		if(curr == 0){
			curr = index;
		}
		Edge temp = prev;
		edges[curr-1] = e;
		edges[curr] = temp; 
		//iter_swap(ptr, prv);
		//ptr--;
		//prv--;
		curr--;
		prev = edges[index-1];
		e = edges[index];
		//e = *ptr;
	}
	*/
	std::vector<Edge>::iterator p;
	p = edges.begin() + index;
	std::sort(edges.begin(), p, pathSort);
	/*
	Edge e = edges[index];
	Edge prev = edges[index-1];
	Edge temp;
	while(e.curr_x < prev.curr_x){
		temp = e;
		edges[curr] = prev;
		edges[curr-1] = temp;
		curr--;
		if(curr == 0){
			break;
		}
		e = edges[curr];
		prev = edges[curr-1];
	}
	*/
	//index = curr;
	
}
static inline bool pathSort(const Edge& e1, const Edge& e2){
		return e1.curr_x < e2.curr_x;
}
static inline bool edgeSort(const Edge& e1, const Edge& e2){
	if(e1.top_y != e2.top_y){
		return e1.top_y < e2.top_y;
	}
	else if(e1.curr_x != e2.curr_x){
		return e1.curr_x < e2.curr_x;
	}
	else{
		return e1.slope < e2.slope;
	}
}





static inline void clipper(const GPoint& p0, const GPoint& p1, std::vector <Edge>& edges, const GRect& clip){
	Edge e1;
	Edge e2;
	Edge e3;
	GPoint p2;
	GPoint  p3;                          
	GPoint q0;
	GPoint q1;
	int w_mult = 1;
	if(GRoundToInt(p0.y()) == GRoundToInt(p1.y())){
		return;
	}

	q0 = p0;
	q1 = p1;
		if(q1.y() < q0.y()){
		GPoint temp1 = q0;
			q0 = q1;
			q1 = temp1;
			w_mult = w_mult * -1;
		}
		//w_mult = w_mult * -1;
			if(q1.y() <= clip.top() || q0.y() >= clip.bottom()) {
		return;
	}
		if(q0.y() < clip.top()){
			double x = ((q1.x()-q0.x())/(q1.y()-q0.y()))*(clip.top() - q0.y());
			q0.set(q0.x()+x, clip.top());
			}
		if(q1.y() > clip.bottom()){
			double x = ((q1.x()-q0.x())/(q1.y()-q0.y()))*(clip.bottom() - q1.y());
			q1.set(q1.x()+x, clip.bottom());
		}
		
		
		
		//left and right
		if(q0.x() > q1.x()){
		GPoint temp2 = q1;
			q1 = q0;
			q0 = temp2;
			w_mult = w_mult*-1;
		}
		//w_mult = w_mult*-1;
		if(q1.x() <= clip.left()){
			q0.set(clip.left(), q0.y());
			q1.set(clip.left(), q1.y());
			if(makeEdge(e1, q0, q1, w_mult)){
				edges.push_back(e1);
			}
			return;
		}
		if(q0.x() >= clip.right()){
			q0.set(clip.right(), q0.y());
			q1.set(clip.right(), q1.y());
			if(makeEdge(e1, q0, q1, w_mult)){
			edges.push_back(e1);
			}
			return;
		}
		if(q0.x() < clip.left()){
			p2.set(clip.left(), q0.y());
			double y = q0.y() + ((q1.y()-q0.y())/(q1.x()-q0.x()))*(clip.left() - q0.x());
			q0.set(clip.left(), y);
			//w_mult = w_mult * -1;
			if(makeEdge(e2, q0, p2, w_mult)){
			edges.push_back(e2);
			}
			
		}
		if(q1.x() > clip.right()){
			p3.set(clip.right(), q1.y());
			double y = q1.y() + ((q1.y()-q0.y())/(q1.x()-q0.x()))*(clip.width() - q1.x());
			q1.set(clip.right(), y);
			
			if(makeEdge(e3, q1, p3, w_mult)){
			edges.push_back(e3); 	
			}
		} 	
			if(makeEdge(e1, q0, q1, w_mult)){
			edges.push_back(e1);
			}
			return;
		
		
		
	
	}
		
	
	


static inline bool makeEdge(Edge& e, const GPoint& p0, const GPoint& p1, int w_mult){
//check floating points
	int top_y;
	int bottom_y;
	if(p0.y() == p1.y() && p0.x() == p1.x()){
		return false;
	}
	
	double slope = 0;
	double curr_x = 0;
	if((p1.y()-p0.y()) == 0){
		return false;
	}
	if(GRoundToInt(p0.y()) < GRoundToInt(p1.y())){
		w_mult*=-1;
		top_y = GRoundToInt(p0.y());
		bottom_y = GRoundToInt(p1.y());	
		curr_x = p0.x();
		slope = (p0.x()-p1.x())/(p0.y()-p1.y());
	}
	else if(GRoundToInt(p1.y()) < GRoundToInt(p0.y())){
		
		top_y = GRoundToInt(p1.y());
		bottom_y = GRoundToInt(p0.y());
		curr_x = p1.x();
		slope = (p1.x()-p0.x())/(p1.y()-p0.y());
	}
	else{
	return false;	
	}
	e.top_y = top_y;
	e.bottom_y = bottom_y;
	e.slope = slope;
	e.curr_x = curr_x;
	e.w = w_mult;
	return true;
}
GRect rectFixer(GRect rect){
	Layer last;
	last = lastLayer(layerStack);
	int top = GRoundToInt(rect.top());
	int left = GRoundToInt(rect.left());
	int right = GRoundToInt(rect.right());
	int bottom = GRoundToInt(rect.bottom());
	if(GRoundToInt(rect.top()) < 0 || GRoundToInt(rect.left()) < 0 || GRoundToInt(rect.right()) > last.bounds.width() || GRoundToInt(rect.bottom()) > last.bounds.height()){
		
		if(GRoundToInt(rect.top()) < 0){
			top = 0;
		}
		if(GRoundToInt(rect.left()) < 0){
			left = 0;
		}
		if(GRoundToInt(rect.right()) > last.bounds.width()){
			right = last.bounds.width();
		}
		if(GRoundToInt(rect.bottom()) > last.bounds.height()){
			bottom = last.bounds.height(); 
		}
		}
		rect.setLTRB(left, top, right, bottom);
		GRect rectFinal = rect;
		return rectFinal;
	}



private:
    const GBitmap fDevice;
	std::vector<Layer> layerStack;
	
static inline GPixel pixelConverter(const GPaint paint){
		GColor color = paint.getColor();
		GColor pinned = color.pinToUnit();
		float a = pinned.fA* 255;
		int a2 = GRoundToInt(a);
		float r = pinned.fR*pinned.fA;
		r = r*255;
		int r2 = GRoundToInt(r);
		float g = pinned.fG*pinned.fA;
		g = g*255;
		int g2 = GRoundToInt(g);
		float b = pinned.fB*pinned.fA;
		b = b*255;
		int b2 = GRoundToInt(b);
		return GPixel_PackARGB(a2, r2, g2, b2);
	} 
static inline GPixel pixelBlender(const GPixel src, const GPixel dest, const GBlendMode blendMode){
			int a = GPixel_GetA(src);
			int r = GPixel_GetR(src);
			int g = GPixel_GetG(src);
			int b = GPixel_GetB(src);
			int a2 = GPixel_GetA(dest);
			int r2 = GPixel_GetR(dest);
			int g2 = GPixel_GetG(dest);
			int b2 = GPixel_GetB(dest);
		if(blendMode == GBlendMode::kClear){ //kClear
			GPixel pixel = GPixel_PackARGB(0, 0, 0, 0);
			return pixel;
			}
		else if(blendMode == GBlendMode::kSrc){ //kSrc
			return src;
		}
		else if(blendMode == GBlendMode::kDst){ //kDst
			return dest;
		}
		else if(blendMode == GBlendMode::kSrcOver){ //kSrcOver
			int aBlend = (255-a);
			aBlend = aBlend*a2;
			aBlend += 127;
			aBlend = aBlend/255;
			int aBlend2 = a + aBlend;
			int rBlend = (255-a);
			rBlend = rBlend *r2;
			rBlend += 127;
			rBlend = rBlend/255;
			int rBlend2 = r + rBlend;
			int gBlend = (255-a);
			gBlend = gBlend*g2;
			gBlend += 127;
			gBlend = gBlend/255;
			int gBlend2= g + gBlend;
			int bBlend = (255-a);
			bBlend = bBlend*b2;
			bBlend += 127;
			bBlend = bBlend/255;
			int bBlend2 = b + bBlend;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kDstOver){ //kDstOver
		//!< [Da + Sa * (1 - Da), Dc + Sc * (1 - Da)]
			int aBlend = (255-a2);
			aBlend = aBlend*a;
			aBlend += 127;
			aBlend = aBlend/255;
			int aBlend2 = a2 + aBlend;
			int rBlend = (255-a2);
			rBlend = rBlend *r;
			rBlend += 127;
			rBlend = rBlend/255;
			int rBlend2 = r2 + rBlend;
			int gBlend = (255-a2);
			gBlend = gBlend*g;
			gBlend += 127;
			gBlend = gBlend/255;
			int gBlend2 = g2 + gBlend;
			int bBlend = (255-a2);
			bBlend = bBlend*b;
			bBlend += 127;
			bBlend = bBlend/255;
			int bBlend2 = b2 + bBlend;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kSrcIn){  //kSrcIn //!< [Sa * Da, Sc * Da] a*a2, c*a2
			//a
			int aBlend = a*a2;
			aBlend += 127;
			int aBlend2 = aBlend/255;
			//r
			int rBlend = r*a2;
			rBlend += 127;
			int rBlend2 = rBlend/255;
			//g
			int gBlend = g*a2;
			gBlend += 127;
			int gBlend2 = gBlend/255;
			//b
			int bBlend = b*a2;
			bBlend += 127;
			int bBlend2 = bBlend/255;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kDstIn){ //kDstIn //!< [Da * Sa, Dc * Sa] a2 * a, c2*a
				//a
			int aBlend = a2*a;
			aBlend += 127;
			int aBlend2 = aBlend/255;
			//r
			int rBlend = r2*a;
			rBlend += 127;
			int rBlend2 = rBlend/255;
			//g
			int gBlend = g2*a;
			gBlend += 127;
			int gBlend2 = gBlend/255;
			//b
			int bBlend = b2*a;
			bBlend += 127;
			int bBlend2 = bBlend/255;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kSrcOut){ //kSrcOut  //!< [Sa * (1 - Da), Sc * (1 - Da)] //a*(255-a2)/255, c*(255-a2)/255
			int aBlend = (255-a2);
			aBlend = aBlend*a;
			aBlend += 127;
			int aBlend2 = aBlend/255;
			int rBlend = (255-a2);
			rBlend = rBlend *r;
			rBlend += 127;
			int rBlend2 = rBlend/255;
			int gBlend = (255-a2);
			gBlend = gBlend*g;
			gBlend += 127;
			int gBlend2 = gBlend/255;
			int bBlend = (255-a2);
			bBlend = bBlend*b;
			bBlend += 127;
			int bBlend2 = bBlend/255;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
		return blend;
		}
		else if(blendMode == GBlendMode::kDstOut){ //kDstOut  //!< [Da * (1 - Sa), Dc * (1 - Sa)] //a2*(255-a)/255, c2 * (255-a)/255
			int aBlend = (255-a);
			aBlend = aBlend*a2;
			aBlend += 127;
			int aBlend2 = aBlend/255;
			int rBlend = (255-a);
			rBlend = rBlend *r2;
			rBlend += 127;
			int rBlend2 = rBlend/255;
			int gBlend = (255-a);
			gBlend = gBlend*g2;
			gBlend += 127;
			int gBlend2 = gBlend/255;
			int bBlend = (255-a);
			bBlend = bBlend*b2;
			bBlend += 127;
			int bBlend2 = bBlend/255;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kSrcATop){ //kSrcATop //!< [Da, Sc * Da + Dc * (1 - Sa)] // a2, ((c*a2) + (c2*(255-a)))/255
			//a
			int aBlend2 = a2;
			//r
			int rBlendF1 = r*a2;
			rBlendF1 += 127;
			rBlendF1 = rBlendF1/255;
			int rBlendF2 = 255-a;
			rBlendF2 = r2*rBlendF2;
			rBlendF2 += 127;
			rBlendF2 = rBlendF2/255;
			int rBlend2 = rBlendF1 + rBlendF2;
			//g
			int gBlendF1 = g*a2;
			gBlendF1 += 127;
			gBlendF1 = gBlendF1/255;
			int gBlendF2 = 255-a;
			gBlendF2 = g2*gBlendF2;
			gBlendF2 += 127;
			gBlendF2 = gBlendF2/255;
			int gBlend2 = gBlendF1 + gBlendF2;
			//b
			int bBlendF1 = b*a2;
			bBlendF1 += 127;
			bBlendF1 = bBlendF1/255;
			int bBlendF2 = 255-a;
			bBlendF2 = b2*bBlendF2;
			bBlendF2 += 127;
			bBlendF2 = bBlendF2/255;
			int bBlend2 = bBlendF1 + bBlendF2;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
			}
		else if(blendMode == GBlendMode::kDstATop){ //kDstATop//!< [Sa, Dc * Sa + Sc * (1 - Da)] //a, ((c2*a) + (c*(255-a2)))/255
			int aBlend2 = a;
			//r
			int rBlendF1 = r2*a;
			rBlendF1 += 127;
			rBlendF1 = rBlendF1/255;
			int rBlendF2 = 255-a2;
			rBlendF2 = r*rBlendF2;
			rBlendF2 += 127;
			rBlendF2 = rBlendF2/255;
			int rBlend2 = rBlendF1 + rBlendF2;
			//g
			int gBlendF1 = g2*a;
			gBlendF1 += 127;
			gBlendF1 = gBlendF1/255;
			int gBlendF2 = 255-a2;
			gBlendF2 = g*gBlendF2;
			gBlendF2 += 127;
			gBlendF2 = gBlendF2/255;
			int gBlend2 = gBlendF1 + gBlendF2;
			//b
			int bBlendF1 = b2*a;
			bBlendF1 += 127;
			bBlendF1 = bBlendF1/255;
			int bBlendF2 = 255-a2;
			bBlendF2 = b*bBlendF2;
			bBlendF2 += 127;
			bBlendF2 = bBlendF2/255;
			int bBlend2 = bBlendF1 + bBlendF2;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kXor){ //kXor //!< [Sa + Da - 2 * Sa * Da, Sc * (1 - Da) + Dc * (1 - Sa)] // ((a + a2) - (510*a*a2))/255, (c*(255-a2) + (c2*(255-a)))/255 
			//a
			int aBlendF1 = a+a2;
			int aBlendF2 = 510*a;
			aBlendF2 += 127;
			aBlendF2 = aBlendF2/255;
			aBlendF2 = a2*aBlendF2;
			aBlendF2 +=127;
			aBlendF2 = aBlendF2/255;
			int aBlend2 = aBlendF1 - aBlendF2;
			//r
			int rBlendF1 = 255-a2;
			rBlendF1 = r*rBlendF1;
			rBlendF1 += 127;
			rBlendF1 = rBlendF1/255;
			int rBlendF2 = 255-a;
			rBlendF2 = r2*rBlendF2;
			rBlendF2 += 127;
			rBlendF2 = rBlendF2/255;
			int rBlend2 = rBlendF1 + rBlendF2;
			//g
			int gBlendF1 = 255-a2;
			gBlendF1 = g*gBlendF1;
			gBlendF1 += 127;
			gBlendF1 = gBlendF1/255;
			int gBlendF2 = 255-a;
			gBlendF2 =g2*gBlendF2;
			gBlendF2 += 127;
			gBlendF2 = gBlendF2/255;
			int gBlend2 = gBlendF1 + gBlendF2;
			//b 
			int bBlendF1 = 255-a2;
			bBlendF1 = b*bBlendF1;
			bBlendF1 += 127;
			bBlendF1 = bBlendF1/255;
			int bBlendF2 = 255-a;
			bBlendF2 = b2*bBlendF2;
			bBlendF2 += 127;
			bBlendF2 = bBlendF2/255;
			int bBlend2 = bBlendF1 + bBlendF2;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}else{
		return GPixel_PackARGB(0, 0, 0, 0);
	}
}

Layer lastLayer(std::vector<Layer> layers){
	//GRect canvas;
	Layer last;
	Layer current;
	for(std::vector<Layer>::reverse_iterator i = layerStack.rbegin(); i != layerStack.rend(); ++i){
			current = *i;
			if(current.isLayer){
				last = current;
				//canvas.setXYWH(0, 0, last.bounds.width(), last.bounds.height());
				break;
		//GBitmap(int w, int h, size_t rb, GPixel* pixels, bool isOpaque)
		}
	}
	return last;
}

/*void layerRect(std::vector<Layer> layers, GRect& rect2){
	//for(vector<Layer>::reverse::iterator i = layerStack.rbegin(); i != layerStack.rend(); ++i){
	//if(i.isLayer()){
		Layer last = lastLayer(layers);
		GRect canvas;
		canvas.setXYWH(0, 0, last.bounds.width(), last.bounds.height());
		//rect2 = rectFixer(rect2);
		 if(rect2.intersects(canvas)){
		 //setLTRB(T l, T t, T r, T b)
		 rect2 = rectFixer(rect2);
		 	return;
		 }else{
		 	rect2.setLTRB(0, 0, 0, 0);	
		 	return; 
		 }
		
		//break;
		//GBitmap(int w, int h, size_t rb, GPixel* pixels, bool isOpaque)
	//}
	//}
}
*/

//end of EMPTYCANVAS
};









std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& device) {
    if (!device.pixels()) {
        return nullptr;
    }
    return std::unique_ptr<GCanvas>(new EmptyCanvas(device));
}
	

  /**
     *  Modifies the CTM by preconcatenating the specified matrix with the CTM. The canvas
     *  is constructed with an identity CTM.
     *
     *  CTM' = CTM * matrix
     */

/*

  GMatrix& preConcat(const GMatrix& primo) {
        this->setConcat(*this, primo);
        return *this;
    }
*/

 class MatrixShader : public GShader {
 	
    public:
        MatrixShader(const GBitmap& device, const GMatrix& localM, GShader::TileMode tm) : sDevice(device), lM(localM), stm(tm) {}
    
	private: 
		const GBitmap sDevice;
		const GMatrix lM;    
		GMatrix current;
		GMatrix inverse;
		GShader::TileMode stm;
		
     bool isOpaque() override{
		return this->sDevice.isOpaque();
	}
     bool setContext(const GMatrix& ctm) override{
		//GMatrix temp;
		//GMatrix localMatrix;
		 //lInv.invert(&localMatrix);
		//temp.setConcat(ctm, localMatrix);
		if(!ctm.invert(&current)){
			return false;
		}
		//9-28 inverted lInv
		lM.invert(&inverse);
		current.postConcat(inverse);
		current.postScale(1.0f/(sDevice.width()), 1.0f/(sDevice.height()));
		return true;
	}
		
        void shadeRow(int x, int y, int count, GPixel row[]) override{
			
for(int i = 0; i < count; i++){
			GPoint localPoint = current.mapXY(x+i+0.5, y+0.5);
			float new_x = localPoint.x();
	 		float new_y = localPoint.y();
			if(stm == GShader::TileMode::kClamp){
			 	new_x = std::fmin(std::fmax(localPoint.fX, 0), 0.99999);
				new_y = std::fmin(std::fmax(localPoint.fY, 0),  0.99999);
			 }
			 if(stm == GShader::TileMode::kRepeat){
			 	if(new_x < 0 || new_x > 1){
			 	new_x = localPoint.fX - floor(localPoint.fX);
				}
				if(new_y < 0 || new_y > 1){ 
				new_y = localPoint.fY - floor(localPoint.fY);
				}
			 }
			 if(stm == GShader::TileMode::kMirror){
			 if(new_x < 0 || new_x > 1){
			 	//new_x = localPoint.fX;
			 	new_x = new_x * 0.5;
				new_x = new_x -floor(new_x);
				if(new_x > 0.5){
					new_x = 1-new_x;
				}
				new_x = new_x*2;
				}
				if(new_y < 0 || new_y > 1){
				//new_y = localPoint.fY;
				//new_y = new_y *2;
				new_y = new_y * 0.5;
				new_y = new_y -floor(new_y);
				if(new_y > 0.5){
					new_y = 1-new_y;
				}
				new_y = new_y *2;
				}
			 }
			 new_x *= this->sDevice.width();
			 new_y *= this->sDevice.height(); 
			// current.postScale(localPoint.x(), localPoint.y());
			 localPoint = GPoint::Make(new_x, new_y);
			/*if(localPoint.x() >= this->sDevice.width()){
					localPoint.fX = this->sDevice.width()-1;
				}
				if(localPoint.y() >= this->sDevice.height()){
					localPoint.fY = this->sDevice.height() - 1;
				}
				if(localPoint.x() < 0){
					localPoint.fX = 0;
				}
				if(localPoint.y() < 0){
					localPoint.fY = 0;
				}*/
			//	current.postScale(this->sDevice.width(), this->sDevice.height());
				GPixel* pixel = this->sDevice.getAddr(new_x, new_y);
				row[i] = *pixel;
				//localPoint.fX += current[GMatrix::SX];
				//localPoint.fY += current[GMatrix::KY];
				
			}
		}
    };
	

 /*   std::unique_ptr<GShader> GCreateBitmapShader(const GBitmap& device, const GMatrix& localInv){
		if(!device.pixels()){
			return nullptr;
		}
		return std::unique_ptr<GShader>(new MatrixShader(device, localInv));
		//device w/localInv
	}
	
/*	class CanvasFilter : public GFilter {
		public:
			CanvasFilter(GBlendMode mode, const GColor& src) : sMode(mode), cSrc(src) {}
			
			private: 
				const GBlendMode sMode;
				const GColor cSrc;    
				GPixel currentOut[];
				GPixel currentIn[];
				int curr;
				
				
	        
	            virtual bool preservesAlpha() override{
					for(int i = 0; i < curr; i++){
						if(currentOut[i] == currentIn[i]){
							continue;
						}
						return false;
					}
					return true;	
				}
				
/*
 *  Apply the specified blendmode and src color to the filter colors. This is similar to
 *  drawing with the src color onto the input colors.
 *
 *  e.g. in filter(output, input, count)
 *      loop(count)
 *         output[i] = blend(src, input[i])
 */
 //static inline GPixel pixelBlender(const GPixel src, const GPixel dest, const GBlendMode blendMode)
	/*            virtual void filter(GPixel output[], const GPixel input[], int count) override{
					for(int i = 0; i < count; i++){
						GPixel cPixel = pixelConverter(cSrc);
						output[i] = pixelBlender(cPixel, input[i], sMode);
					}
				}
static inline GPixel pixelConverter(const GColor color){
		//GColor color = paint.getColor();
		GColor pinned = color.pinToUnit();
		float a = pinned.fA* 255;
		int a2 = GRoundToInt(a);
		float r = pinned.fR*pinned.fA;
		r = r*255;
		int r2 = GRoundToInt(r);
		float g = pinned.fG*pinned.fA;
		g = g*255;
		int g2 = GRoundToInt(g);
		float b = pinned.fB*pinned.fA;
		b = b*255;
		int b2 = GRoundToInt(b);
		return GPixel_PackARGB(a2, r2, g2, b2);
} 
static inline GPixel pixelBlender(const GPixel src, const GPixel dest, const GBlendMode blendMode){
			int a = GPixel_GetA(src);
			int r = GPixel_GetR(src);
			int g = GPixel_GetG(src);
			int b = GPixel_GetB(src);
			int a2 = GPixel_GetA(dest);
			int r2 = GPixel_GetR(dest);
			int g2 = GPixel_GetG(dest);
			int b2 = GPixel_GetB(dest);
		if(blendMode == GBlendMode::kClear){ //kClear
			GPixel pixel = GPixel_PackARGB(0, 0, 0, 0);
			return pixel;
			}
		else if(blendMode == GBlendMode::kSrc){ //kSrc
			return src;
		}
		else if(blendMode == GBlendMode::kDst){ //kDst
			return dest;
		}
		else if(blendMode == GBlendMode::kSrcOver){ //kSrcOver
			int aBlend = (255-a);
			aBlend = aBlend*a2;
			aBlend += 127;
			aBlend = aBlend/255;
			int aBlend2 = a + aBlend;
			int rBlend = (255-a);
			rBlend = rBlend *r2;
			rBlend += 127;
			rBlend = rBlend/255;
			int rBlend2 = r + rBlend;
			int gBlend = (255-a);
			gBlend = gBlend*g2;
			gBlend += 127;
			gBlend = gBlend/255;
			int gBlend2= g + gBlend;
			int bBlend = (255-a);
			bBlend = bBlend*b2;
			bBlend += 127;
			bBlend = bBlend/255;
			int bBlend2 = b + bBlend;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kDstOver){ //kDstOver
		//!< [Da + Sa * (1 - Da), Dc + Sc * (1 - Da)]
			int aBlend = (255-a2);
			aBlend = aBlend*a;
			aBlend += 127;
			aBlend = aBlend/255;
			int aBlend2 = a2 + aBlend;
			int rBlend = (255-a2);
			rBlend = rBlend *r;
			rBlend += 127;
			rBlend = rBlend/255;
			int rBlend2 = r2 + rBlend;
			int gBlend = (255-a2);
			gBlend = gBlend*g;
			gBlend += 127;
			gBlend = gBlend/255;
			int gBlend2 = g2 + gBlend;
			int bBlend = (255-a2);
			bBlend = bBlend*b;
			bBlend += 127;
			bBlend = bBlend/255;
			int bBlend2 = b2 + bBlend;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kSrcIn){  //kSrcIn //!< [Sa * Da, Sc * Da] a*a2, c*a2
			//a
			int aBlend = a*a2;
			aBlend += 127;
			int aBlend2 = aBlend/255;
			//r
			int rBlend = r*a2;
			rBlend += 127;
			int rBlend2 = rBlend/255;
			//g
			int gBlend = g*a2;
			gBlend += 127;
			int gBlend2 = gBlend/255;
			//b
			int bBlend = b*a2;
			bBlend += 127;
			int bBlend2 = bBlend/255;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kDstIn){ //kDstIn //!< [Da * Sa, Dc * Sa] a2 * a, c2*a
				//a
			int aBlend = a2*a;
			aBlend += 127;
			int aBlend2 = aBlend/255;
			//r
			int rBlend = r2*a;
			rBlend += 127;
			int rBlend2 = rBlend/255;
			//g
			int gBlend = g2*a;
			gBlend += 127;
			int gBlend2 = gBlend/255;
			//b
			int bBlend = b2*a;
			bBlend += 127;
			int bBlend2 = bBlend/255;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kSrcOut){ //kSrcOut  //!< [Sa * (1 - Da), Sc * (1 - Da)] //a*(255-a2)/255, c*(255-a2)/255
			int aBlend = (255-a2);
			aBlend = aBlend*a;
			aBlend += 127;
			int aBlend2 = aBlend/255;
			int rBlend = (255-a2);
			rBlend = rBlend *r;
			rBlend += 127;
			int rBlend2 = rBlend/255;
			int gBlend = (255-a2);
			gBlend = gBlend*g;
			gBlend += 127;
			int gBlend2 = gBlend/255;
			int bBlend = (255-a2);
			bBlend = bBlend*b;
			bBlend += 127;
			int bBlend2 = bBlend/255;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
		return blend;
		}
		else if(blendMode == GBlendMode::kDstOut){ //kDstOut  //!< [Da * (1 - Sa), Dc * (1 - Sa)] //a2*(255-a)/255, c2 * (255-a)/255
			int aBlend = (255-a);
			aBlend = aBlend*a2;
			aBlend += 127;
			int aBlend2 = aBlend/255;
			int rBlend = (255-a);
			rBlend = rBlend *r2;
			rBlend += 127;
			int rBlend2 = rBlend/255;
			int gBlend = (255-a);
			gBlend = gBlend*g2;
			gBlend += 127;
			int gBlend2 = gBlend/255;
			int bBlend = (255-a);
			bBlend = bBlend*b2;
			bBlend += 127;
			int bBlend2 = bBlend/255;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kSrcATop){ //kSrcATop //!< [Da, Sc * Da + Dc * (1 - Sa)] // a2, ((c*a2) + (c2*(255-a)))/255
			//a
			int aBlend2 = a2;
			//r
			int rBlendF1 = r*a2;
			rBlendF1 += 127;
			rBlendF1 = rBlendF1/255;
			int rBlendF2 = 255-a;
			rBlendF2 = r2*rBlendF2;
			rBlendF2 += 127;
			rBlendF2 = rBlendF2/255;
			int rBlend2 = rBlendF1 + rBlendF2;
			//g
			int gBlendF1 = g*a2;
			gBlendF1 += 127;
			gBlendF1 = gBlendF1/255;
			int gBlendF2 = 255-a;
			gBlendF2 = g2*gBlendF2;
			gBlendF2 += 127;
			gBlendF2 = gBlendF2/255;
			int gBlend2 = gBlendF1 + gBlendF2;
			//b
			int bBlendF1 = b*a2;
			bBlendF1 += 127;
			bBlendF1 = bBlendF1/255;
			int bBlendF2 = 255-a;
			bBlendF2 = b2*bBlendF2;
			bBlendF2 += 127;
			bBlendF2 = bBlendF2/255;
			int bBlend2 = bBlendF1 + bBlendF2;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
			}
		else if(blendMode == GBlendMode::kDstATop){ //kDstATop//!< [Sa, Dc * Sa + Sc * (1 - Da)] //a, ((c2*a) + (c*(255-a2)))/255
			int aBlend2 = a;
			//r
			int rBlendF1 = r2*a;
			rBlendF1 += 127;
			rBlendF1 = rBlendF1/255;
			int rBlendF2 = 255-a2;
			rBlendF2 = r*rBlendF2;
			rBlendF2 += 127;
			rBlendF2 = rBlendF2/255;
			int rBlend2 = rBlendF1 + rBlendF2;
			//g
			int gBlendF1 = g2*a;
			gBlendF1 += 127;
			gBlendF1 = gBlendF1/255;
			int gBlendF2 = 255-a2;
			gBlendF2 = g*gBlendF2;
			gBlendF2 += 127;
			gBlendF2 = gBlendF2/255;
			int gBlend2 = gBlendF1 + gBlendF2;
			//b
			int bBlendF1 = b2*a;
			bBlendF1 += 127;
			bBlendF1 = bBlendF1/255;
			int bBlendF2 = 255-a2;
			bBlendF2 = b*bBlendF2;
			bBlendF2 += 127;
			bBlendF2 = bBlendF2/255;
			int bBlend2 = bBlendF1 + bBlendF2;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}
		else if(blendMode == GBlendMode::kXor){ //kXor //!< [Sa + Da - 2 * Sa * Da, Sc * (1 - Da) + Dc * (1 - Sa)] // ((a + a2) - (510*a*a2))/255, (c*(255-a2) + (c2*(255-a)))/255 
			//a
			int aBlendF1 = a+a2;
			int aBlendF2 = 510*a;
			aBlendF2 += 127;
			aBlendF2 = aBlendF2/255;
			aBlendF2 = a2*aBlendF2;
			aBlendF2 +=127;
			aBlendF2 = aBlendF2/255;
			int aBlend2 = aBlendF1 - aBlendF2;
			//r
			int rBlendF1 = 255-a2;
			rBlendF1 = r*rBlendF1;
			rBlendF1 += 127;
			rBlendF1 = rBlendF1/255;
			int rBlendF2 = 255-a;
			rBlendF2 = r2*rBlendF2;
			rBlendF2 += 127;
			rBlendF2 = rBlendF2/255;
			int rBlend2 = rBlendF1 + rBlendF2;
			//g
			int gBlendF1 = 255-a2;
			gBlendF1 = g*gBlendF1;
			gBlendF1 += 127;
			gBlendF1 = gBlendF1/255;
			int gBlendF2 = 255-a;
			gBlendF2 =g2*gBlendF2;
			gBlendF2 += 127;
			gBlendF2 = gBlendF2/255;
			int gBlend2 = gBlendF1 + gBlendF2;
			//b 
			int bBlendF1 = 255-a2;
			bBlendF1 = b*bBlendF1;
			bBlendF1 += 127;
			bBlendF1 = bBlendF1/255;
			int bBlendF2 = 255-a;
			bBlendF2 = b2*bBlendF2;
			bBlendF2 += 127;
			bBlendF2 = bBlendF2/255;
			int bBlend2 = bBlendF1 + bBlendF2;
			GPixel blend = GPixel_PackARGB(aBlend2, rBlend2, gBlend2, bBlend2);
			return blend;
		}else{
		return GPixel_PackARGB(0, 0, 0, 0);
	}
}

};
*/

class ColorGradient : public GShader {
	public:
		ColorGradient(GPoint p0, GPoint p1, const GColor colors[], int count, GShader::TileMode tm) : po0(p0), po1(p1), sCount(count), stm(tm){
		sColors = new GColor[count]; 
		memcpy(sColors, colors, count*sizeof(GColor));
		}
	private:
		GPoint po0;
		GPoint po1;
		GColor* sColors;
		int sCount;	
		GMatrix current;
		GShader::TileMode stm;
		
	 bool isOpaque() override{
	 	for(int i = 0; i < sCount; i++){
		GPixel curr = pixelConverter(sColors[i]);
		if(GPixel_GetA(curr) == 1){
			continue;
		}
		return false;
		}
		return true;
	}
     bool setContext(const GMatrix& ctm) override{
	 GPoint pZero;
	 GPoint pOne;
	double dx = po1.x() - po0.x();
	double dy = po1.y() - po0.y();
	double dist = sqrt((dx*dx) + (dy*dy));
	 GMatrix transform;
	 transform.set6(dx, -dy, po0.x(), dy, dx, po0.y());
		GMatrix temp;
		GMatrix localMatrix;
		localMatrix = transform;
		GMatrix lInv;
		localMatrix.invert(&lInv);
		// lInv.invert(&localMatrix);
		//temp.setConcat(ctm, localMatrix);
		if(!ctm.invert(&current)){
			return false;
		}
		current.postConcat(lInv);
		return true;
	}
		
        void shadeRow(int x, int y, int count, GPixel row[]) override{

			 for(int i = 0; i < count; i++){
			 	GPoint localPoint = current.mapXY(x + 0.5 + i, y + 0.5);
				//GPoint localPoint = current.mapXY(x+i+0.5, y+0.5);
			float new_x = localPoint.x();
	 		float new_y = localPoint.y();
			if(stm == GShader::TileMode::kClamp){
			 	new_x = std::fmin(std::fmax(localPoint.fX, 0), 0.99999);
				new_y = std::fmin(std::fmax(localPoint.fY, 0),  0.99999);
			 }
			 if(stm == GShader::TileMode::kRepeat){
			 	if(new_x < 0 || new_x > 1){
			 	new_x = localPoint.fX - floor(localPoint.fX);
				}
				if(new_y < 0 || new_y > 1){ 
				new_y = localPoint.fY - floor(localPoint.fY);
				}
			 }
			 if(stm == GShader::TileMode::kMirror){
			 if(new_x < 0 || new_x > 1){
			 	//new_x = localPoint.fX;
			 	new_x = new_x * 0.5;
				new_x = new_x -floor(new_x);
				if(new_x > 0.5){
					new_x = 1-new_x;
				}
				new_x = new_x*2;
				}
				if(new_y < 0 || new_y > 1){
				//new_y = localPoint.fY;
				//new_y = new_y *2;
				new_y = new_y * 0.5;
				new_y = new_y -floor(new_y);
				if(new_y > 0.5){
					new_y = 1-new_y;
				}
				new_y = new_y *2;
				}
			 }
		//	 new_x *= this->sDevice.width();
		//	 new_y *= this->sDevice.height(); 
			// current.postScale(localPoint.x(), localPoint.y());
			localPoint = GPoint::Make(new_x, new_y);
				//old begin
				float t = localPoint.x();
				//printf("t: %f", t);
				if(t < 0){
					t = 0;
				}
				if(t > 1){
					t = 1;
				}
				if(t == 0){
					row[i] = pixelConverter(sColors[0]);
					continue;
				}
				if(t == 1){
					row[i] = pixelConverter(sColors[sCount-1]);
					continue;
				}
				int index = floor(t*(sCount-1));
				//printf("sCount: %i", sCount);
				//printf("index: %i", index);
				GColor c1 = sColors[index];
				GColor c2 = sColors[index + 1];
				float d1 = sCount -1;
			 	d1 =  1/d1;
				float sDst = sCount -1;
				sDst = t*sDst;
				sDst = floor(sDst);
				sDst = sDst *d1;
				sDst = t - sDst;
			 	//float sDst = t - (floor(t*(sCount-1))*d1);
				t = sDst/d1;
				GColor gradient;
				//color 1 
				float col1 = 1-t;
				float a = (c1.fA * col1) + (c2.fA*(t));
				float r = (c1.fR * col1) + (c2.fR*(t));
				float g = (c1.fG * col1) + (c2.fG*(t));
				float b = (c1.fB * col1) + (c2.fB*(t));
				//float a = 255;
				//float r = 0;
				//float g = 255;
				//float b = 0;
				//GColor gradient;
				//GColor c1 = sColors[0];
				//GColor c2 = sColors[1];
				//float a = (c1.fA * (1-t)) + (c2.fA*(t));
				//float r = (c1.fR * (1-t)) + (c2.fR*(t));
				//float g = (c1.fG * (1-t)) + (c2.fG*(t));
				//float b = (c1.fB * (1-t)) + (c2.fB*(t));
				//float a = 1;
				//float r = 0;
				//float g = 1;
				//float b = 0;
				//gradient.MakeARGB(a, r, g, b);
				gradient.fA = a;
				gradient.fR = r;
				gradient.fG = g;
				gradient.fB = b;
				row[i] = pixelConverter(gradient);
			 }
			 //count = number of colors
			 //multiple colors(>2): always evenly spaced
			 //get local distance from left, divide it by distance between the left and right color
			 //d1 = 1/(colorCount-1) 
			 //t - 2*d1 d1 -> distance between any given two colors
			 //2 -> floor(t*(count-1))
			 
			
		} 
	 
	static inline GPixel pixelConverter(const GColor color){
		//GColor color = paint.getColor();
		GColor pinned = color.pinToUnit();
		float a = pinned.fA* 255;
		int a2 = GRoundToInt(a);
		float r = pinned.fR*pinned.fA;
		r = r*255;
		int r2 = GRoundToInt(r);
		float g = pinned.fG*pinned.fA;
		g = g*255;
		int g2 = GRoundToInt(g);
		float b = pinned.fB*pinned.fA;
		b = b*255;
		int b2 = GRoundToInt(b);
		return GPixel_PackARGB(a2, r2, g2, b2);
	} 	
	
	/*Layer lastLayer(std::vector<Layer> layers){
	//GRect canvas;
	Layer last;
	Layer current;
	for(std::vector<Layer>::reverse_iterator i = layerStack.rbegin(); i != layerStack.rend(); ++i){
			current = *i;
			if(current.isLayer){
				last = current;
				//canvas.setXYWH(0, 0, last.bounds.width(), last.bounds.height());
				break;
		//GBitmap(int w, int h, size_t rb, GPixel* pixels, bool isOpaque)
		}
	}
	return last;
}		
*/		

	
};		
        
 /*    std::unique_ptr<GFilter> GCreateBlendFilter(GBlendMode mode, const GColor& src){
			return std::unique_ptr<GFilter>(new CanvasFilter(mode, src));
		}
	*/	
	 std::unique_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor colors[], int count, GShader::TileMode tm){
	 	//return std::unique_ptr<GShader>(new ColorGradient(p0, p1, colors, count));
		return std::unique_ptr<GShader>(new ColorGradient(p0, p1, colors, count, tm));     
	 }
	 
	 void GDrawSomething_polys(GCanvas* canvas){
	GPoint points[3]; 
	GPoint point1;
	point1.set(25,25);
	GPoint point2;
	point2.set(100, 100);
	GPoint point3;
	point3.set(50, 150);
	points[0] = point1;
	points[1] = point2;
	points[2] = point3;
	canvas->drawConvexPolygon(points, 3, GPaint({1, 0, 1, 0}));
}


 /*class MatrixShader : public GShader {
 	
    public:
        MatrixShader(const GBitmap& device, const GMatrix& localInv, GShader::TileMode tm) : sDevice(device), lInv(localInv), stm(tm) {}
    
	private: 
		const GBitmap sDevice;
		const GMatrix lInv;    
		GShader::TileMode stm;
		GMatrix current;
		
     bool isOpaque() override{
		return this->sDevice.isOpaque();
	}
     bool setContext(const GMatrix& ctm) override{
		//GMatrix temp;
		//GMatrix localMatrix;
		 //lInv.invert(&localMatrix);
		//temp.setConcat(ctm, localMatrix);
		if(!ctm.invert(&current)){
			return false;
		}
		current.postConcat(lInv);
		//if(stm == GShader::TileMode::kClamp){
		current.postScale(1.0f/(sDevice.width()), 1.0f/(sDevice.height()));
		
		//}
		
		return true;
	}
		
        void shadeRow(int x, int y, int count, GPixel row[]) override{
	/**
     *  Given a row of pixels in device space [x, y] ... [x + count - 1, y], return the
     *  corresponding src pixels in row[0...count - 1]. The caller must ensure that row[]
     *  can hold at least [count] entries.
     */
	 
	 

	 //scale up
	 
			//GMatrix fake = current;
			//fake.setScale(this->sDevice.width(), this->sDevice.height());
	/*		for(int i = 0; i < count; i++){
			GPoint localPoint = current.mapXY(x+i+0.5, y+0.5);
			float new_x = localPoint.x();
	 		float new_y = localPoint.y();
			if(stm == GShader::TileMode::kClamp){
			 	new_x = std::fmin(std::fmax(localPoint.fX, 0), 0.99999);
				new_y = std::fmin(std::fmax(localPoint.fY, 0),  0.99999);
			 }
			 if(stm == GShader::TileMode::kRepeat){
			 	if(new_x < 0 || new_x > 1){
			 	new_x = localPoint.fX - floor(localPoint.fX);
				}
				if(new_y < 0 || new_y > 1){ 
				new_y = localPoint.fY - floor(localPoint.fY);
				}
			 }
			 if(stm == GShader::TileMode::kMirror){
			 if(new_x < 0 || new_x > 1){
			 	//new_x = localPoint.fX;
			 	new_x = new_x * 0.5;
				new_x = new_x -floor(new_x);
				if(new_x > 0.5){
					new_x = 1-new_x;
				}
				new_x = new_x*2;
				}
				if(new_y < 0 || new_y > 1){
				//new_y = localPoint.fY;
				//new_y = new_y *2;
				new_y = new_y * 0.5;
				new_y = new_y -floor(new_y);
				if(new_y > 0.5){
					new_y = 1-new_y;
				}
				new_y = new_y *2;
				}
			 }
			 new_x *= this->sDevice.width();
			 new_y *= this->sDevice.height(); 
			// current.postScale(localPoint.x(), localPoint.y());
			 localPoint = GPoint::Make(new_x, new_y);
			/*if(localPoint.x() >= this->sDevice.width()){
					localPoint.fX = this->sDevice.width()-1;
				}
				if(localPoint.y() >= this->sDevice.height()){
					localPoint.fY = this->sDevice.height() - 1;
				}
				if(localPoint.x() < 0){
					localPoint.fX = 0;
				}
				if(localPoint.y() < 0){
					localPoint.fY = 0;
				}*/
			//	current.postScale(this->sDevice.width(), this->sDevice.height());
			//	GPixel* pixel = this->sDevice.getAddr(new_x, new_y)
			//	row[i] = *pixel;
				//localPoint.fX += current[GMatrix::SX];
				//localPoint.fY += current[GMatrix::KY];
				
	//		}
	//	}
    //};
	
	    std::unique_ptr<GShader> GCreateBitmapShader(const GBitmap& device, const GMatrix& localInv, GShader::TileMode tm){
		if(!device.pixels()){
			return nullptr;
		}
		return std::unique_ptr<GShader>(new MatrixShader(device, localInv, tm));
		//device w/localInv
	}

