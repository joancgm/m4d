
typedef struct array_
{ 
  char *name;
  int size;        /* size in d,f,i,c or p terms */
  char type;         /* f-float d-double i-int c-char p-pointers   */
  void *pointer;   /*address of array values */
  struct array_ *next;    /* next array =0 if none */
  struct array_ *previous;  /* previous array =0 if none */
  int alias;    /* number of aliases  set to -1 for aliases */
  struct array_ **aliasarray; /* corresponding array pointers */
} Array;

typedef struct 
{ 	
  void *pointer;   /*address of array values */
  void *previous;  /* previous array =0 if none */
} TArray;


int arraydelete(const char *name);
int arraysize(const char *name);
char arraytype(const char *name);
void *createarray(const char *name, int size, char type, int nynew);
void *find(const char *name);
Array *findarray(const char *name);
void *need(const char *name);
void *smalloc(int i); 
void  *smalloca(int i, char c);   
void *tmalloca(int i, char c);  

