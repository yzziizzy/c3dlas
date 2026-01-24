#ifndef __C3DLAS__internal_array_h__
#define __C3DLAS__internal_array_h__

	
inline static unsigned long ARR_nextPOT(unsigned long in) {
	
	in--;
	in |= in >> 1;
	in |= in >> 2;
	in |= in >> 4;
	in |= in >> 8;
	in |= in >> 16;
	in++;
	
	return in;
}
	
	
#define ARR(type) \
	struct { \
		unsigned long alloc, count; \
		type* data; \
	} 


#define ARR_check(x, sz) \
	do { \
		if((x)->alloc < (sz)) { \
			(x)->alloc = (x)->alloc == 0 ? 8 : ARR_nextPOT(sz); \
			(x)->data = C3DLAS_realloc((x)->data, (x)->alloc * sizeof(*(x)->data)); \
		} \
	} while(0);
	
	
#define ARR_inc(x, out) \
	do { \
		ARR_check(x, (x)->count + 1); \
		(out) = &(x)->data[(x)->count++]; \
	} while(0);
	

// the second argument is the item to push, but varargs to swallow commas in objects without extra parens
#define ARR_push(x, ...) \
	do { \
		ARR_check(x, (x)->count + 1); \
		(x)->data[(x)->count++] = (__VA_ARGS__); \
	} while(0);
	


#define ARR_free(x) \
	do { \
		if((x)->data) { \
			C3DLAS_free((x)->data); \
		} \
		(x)->data = NULL; \
		(x)->count = 0; \
		(x)->alloc = 0; \
	} while(0);


#define ARR_EACH(x, i) for(unsigned long i = 0; i < (x)->count; i++)


#endif // __C3DLAS__internal_array_h__
