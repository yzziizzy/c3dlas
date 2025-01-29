
#include "../c3dlas.h"


int rasterizeLine2d(vec2 pa, vec2 pb, float cellSize, int (*found)(void* userData, s64 x, s64 y), void* userData) {

	float cs = cellSize;
	
	f32 rise = (pb.y - pa.y);
	f32 run = (pb.x - pa.x);
	
	f32 absminy = MIN(pb.y, pa.y);
	f32 absmaxy = MAX(pb.y, pa.y);
	
	if(run == 0) {
		s64 x = floor(pa.x / cs);
		
		s64 miny = floor(absminy / cs);
		s64 maxy = floor(absmaxy / cs);
		
		for(s64 y = miny; y <= maxy; y++) {
			if(found(userData, x, y)) return 1;
		}
	}
	else {
		// y = slope * x + b
		f32 slope = rise / run;
		f32 b = pa.y - slope * pa.x;
		
		// get x limits of the line
		s64 minx = floor(MIN(pa.x, pb.x) / cs);
		s64 maxx = floor(MAX(pa.x, pb.x) / cs);
		
		for(s64 x = minx; x <= maxx; x++) {
			// get y intercepts in this column, for both sides
			f32 xl = x * cs;
			f32 xh = (x + 1) * cs;
			
			f32 y1 = xl * slope + b;
			f32 y2 = xh * slope + b;
			
			f32 yllf = MAX(absminy, MIN(y1, y2));
			f32 yhlf = MIN(absmaxy, MAX(y1, y2));
			
			// get y limits in this column
			s64 miny = floor(MIN(yllf, yhlf) / cs);
			s64 maxy = floor(MAX(yllf, yhlf) / cs);
			
			for(s64 y = miny; y <= maxy; y++) {
				if(found(userData, x, y)) return 1;
			}
		}
	}
	
	return 0;
}


// conservative overestimation; the ends are square, not conformed to the radius
int rasterizeFatLine2d(vec2 pa, vec2 pb, f32 width, float cellSize, int (*found)(void* userData, s64 x, s64 y), void* userData) {

	float cs = cellSize;
	
	f32 rise = (pb.y - pa.y);
	f32 run = (pb.x - pa.x);

	f32 absminy = MIN(pb.y, pa.y) - width;
	f32 absmaxy = MAX(pb.y, pa.y) + width;
	
	if(run == 0) { // line is straight along the y axis
		int x = floor(pa.x / cs);
		
		int miny = floor(absminy / cs);
		int maxy = floor(absmaxy / cs);
		int minx = floor((pa.x - width) / cs);
		int maxx = floor((pa.x + width) / cs);
		
		for(int x = minx; x <= maxx; x++)
		for(int y = miny; y <= maxy; y++) {
			if(found(userData, x, y)) return 1;
		}
	}
	else {
		
		vec2 dir = vNorm(vSub(pb, pa));
		vec2 perp = {dir.y, -dir.x};
		
		vec2 pa1 = vAdd(pa, vScale(perp, width));
		vec2 pa2 = vAdd(pa, vScale(perp, -width));
		vec2 pb1 = vAdd(pb, vScale(perp, width));
		vec2 pb2 = vAdd(pb, vScale(perp, -width));
		
		// y = slope * x + b
		f32 slope = rise / run;
		f32 b1 = pa1.y - slope * pa1.x; // the parallel lines along the top and bottom of the radius
		f32 b2 = pa2.y - slope * pa2.x;
		
		// get x limits of the line
		int minx = floor((MIN(pa.x, pb.x) - width) / cs);
		int maxx = floor((MAX(pa.x, pb.x) + width) / cs);
		
		for(int x = minx; x <= maxx; x++) {
			// get y intercepts in this column, for both sides
			f32 xl = x * cs;
			f32 xh = (x + 1) * cs;
			
			f32 y1a = xl * slope + b1;
			f32 y1b = xl * slope + b2;
			f32 y2a = xh * slope + b1;
			f32 y2b = xh * slope + b2;
			
			f32 yllf = MAX(absminy, MIN(MIN(y1a, y2a), MIN(y1b, y2b)));
			f32 yhlf = MIN(absmaxy, MAX(MAX(y1a, y2a), MAX(y1b, y2b)));
			
			// get y limits in this column
			int miny = floor((MIN(yllf, yhlf)) / cs);
			int maxy = floor((MAX(yllf, yhlf)) / cs);
			
			for(int y = miny; y <= maxy; y++) {
				if(found(userData, x, y)) return 1;
			}	
		}
	}
	
	return 0;
}
	
