/* some algebra */

#define abs(a) ((a)>0?(a):-(a))
#define max(a,b) ((a)>(b)?(a):(b))
#define max3(a,b,c) ((a)>max(b,c)?(a):max(b,c))
#define max4(a,b,c,d) (max(max(a,b),max(c,d)))
#define max8(a,b,c,d,e,f,g,h)(max4(a,b,c,d)>max4(e,f,g,h)?max4(a,b,c,d):max4(e,f,g,h))
#define min(a,b) ((a)<(b)?(a):(b))
#define min4(a,b,c,d) (min(min(a,b),min(c,d)))
