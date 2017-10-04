/* contains p_outgif p_newline p_gifheader p_lzw p_set2char p_writebits */

#include "global.h"
#include "p_plot.h"
/* conver an image to a .gif file (outgif) and routines it uses */
/*-----------------------------------------------------------*/
void p_outgif(char*filename, int*image, int*ipix)
{
  char*rgb; int numpix; 
  FILE *fp=NULL;
  char**defaultdir,*colorlistname,*rgblocal;
  int i,j,ii,ir,ig,ib;
  int *ic;
  
  ic=(int *)tmalloca(ipix[0]*ipix[1],'i');
  /* turn upside down */
  for (i=0;i<ipix[0];i++) for (j=0;j<ipix[1];j++) 
    ic[i+ipix[0]*(ipix[1]-1-j)]=image[i+ipix[0]*j];
  
  /* get rgb colormap if not already read */
  rgb=(char *)find("p.rgbcolor");
  numpix=arraysize("p.rgbcolor")/3;
  if (rgb==0) /* not found - read it */
  {
    printout("normal","looking for file color.map\n");
    fp=fopen("color.map","r");
    if (fp==NULL)
    {
      printout("normal","file color.map not found, try default directory\n");
      defaultdir=(char**)find("p.defaultdir");
      if (defaultdir==0) 
      {
        printout("error p_outgif"," string variable p.defaultdir not found\n no colormap, cannot create gif file, exit\n");
        exitm4d(0);
      }
      for (i=0;i<arraysize("p.defaultdir");i++)
      {
        if (defaultdir[i][strlen(defaultdir[i])-1]=='/')
          colorlistname=setname2(defaultdir[i],"color.map");
        else colorlistname=setname3(defaultdir[i],"/","color.map");  
        fp=fopen(colorlistname,"r");
        if (fp!=NULL) 
        {
          printout("normal","using file %s\n",colorlistname);
          break;
        }
      }
      if (fp==NULL)
      {
        printout("error p_outgif"," no colormap, cannot create gif file, exit\n");
        exitm4d(0);
      }
    }
    /* have file read color.map */
    numpix=readint(fp); 
    p_newline(fp);
    rgb=(char*)createarray("p.rgbcolor",numpix*3,'c',0);
    rgblocal=rgb;
    for (i=0; i< numpix; i++)
    { ii=readint(fp);
      ir=readint(fp);
      ig=readint(fp);
      ib=readint(fp);
      p_newline(fp); 
      *rgblocal= ir; rgblocal++;
      *rgblocal= ig; rgblocal++;
      *rgblocal= ib; rgblocal++;
    }
    fclose(fp);
  } /* end read colormap */
  /* write gif file */
  fp=safefopen(filename,"w");
  p_gifheader(fp,ipix[0], ipix[1], numpix, rgb);   /* gif header             */
  p_lzw(fp,numpix, ipix[0]*ipix[1], ic);         /* lzw compression codes  */
  fprintf(fp,"%c;",(char)0);                    /* end of gif file        */
  fclose(fp);
}

/*-----------------------------------------------------------*/
int p_newline(FILE *fv)
{ 
  int i; char c;
  i=0; c=' ';
  while (c != '\n' && c != EOF) { c=fgetc(fv); i++;}
  return i; 
}
/*-----------------------------------------------------------*/
void p_set2char(int i,char c2[2])
{ 
  c2[1] = i >>8;
  c2[0] = i % 256;
}
/*-----------------------------------------------------------*/
void p_gifheader(FILE *fp,int ixtot,int iytot,int numpix,char *rgb)
{ 
  int i, itot, j;
  char sd[7], id[10];
  
  fprintf(fp,"GIF87a");           /* header */
  p_set2char(ixtot,&sd[0]);
  p_set2char(iytot,&sd[2]);
  j=0; itot=2; while(itot<numpix) {j++; itot *= 2;}
  sd[4] = j + (15 << 4);
  sd[5] = 0;
  sd[6] = 0;
  id[0] = ',';
  for (i=1;i<10;i++) id[i]=0;
  p_set2char(ixtot,&id[5]);
  p_set2char(iytot,&id[7]);
  for (i=0;i<7;i++) fprintf(fp,"%c",sd[i]);  /* logical screen descriptor */
  for (i=0; i<numpix * 3; i++) fprintf(fp,"%c",rgb[i]);  /* colors */
  for (i=numpix*3; i< itot*3; i++) fprintf(fp,"%c",(char)0);
  for (i=0;i<10;i++) fprintf(fp,"%c",id[i]);  /* image descriptor */
} 
/*-----------------------------------------------------------*/

void p_writebits(int iv, int nbits, FILE *fp, int *sizei, int *nextbiti, int outbuf[200])
{
  int size,nextbit;
  int mbits, top, use, i, add;
  size=sizei[0]; nextbit=nextbiti[0];
  if (nbits==0)
  { size=0; nextbit=0;
    for (i=0; i<200; i++) outbuf[i]=0;
  }
  else if (nbits<0)
  { if (size>0)
  { fprintf(fp,"%c", (char)size);
    for (i=0; i<size; i++) fprintf(fp,"%c", (char)outbuf[i]);
  }
  }
  else
    while(nbits>0)
    { 
      if (nextbit==0) 
      { size+=1;
        if (size>200)
        { fprintf(fp,"%c", (char)200);
          for (i=0; i<200; i++) fprintf(fp,"%c", (char)outbuf[i]);
          size=1; 
          for (i=0; i<200; i++) outbuf[i]=0;
        }
      }
      mbits = 8 - nextbit;
      top = iv >> mbits;
      use = iv - (top<<mbits);
      add =  use<<nextbit;
      outbuf[size-1] += add;
      iv=top;
      nextbit +=nbits;
      if (nextbit > 7) nextbit=0;
      nbits -= mbits;
    }
  sizei[0]=size; nextbiti[0]=nextbit;
}
/*-----------------------------------------------------------*/
void p_lzw(FILE *fp,int numpix, int imax, int* ic)
{
  int size=0, nextbit=0, outbuf[200]; /* to pass to writebits */
  int *lzwcode, *prefixcode, *addchar,
  i, j, jtot, itot, nbits, nextcode, itable, clearcode, endcode,
  now, prefix, index, maxbits, tablesize, firstcode, jstart;
  
  /*---------------initialize and set up table space------------*/
  jtot=1; itot=2;  while(itot<numpix) {jtot++; itot *= 2;}
  nbits = jtot+1; /* nbits = 8; */
  itable=0;
  maxbits=12;
  tablesize=2<<maxbits;
  lzwcode=malloc(tablesize * 4);
  prefixcode=malloc(tablesize * 4);
  addchar=malloc(tablesize * 4);
  clearcode=itot; endcode=itot+1; firstcode=itot+2; nextcode=itot+2;
  /*--------------- compress --------------------------*/
  p_writebits(0, 0,fp,&size,&nextbit,outbuf);
  fprintf(fp,"%c",(char)nbits-1);
  p_writebits(clearcode, nbits,fp,&size,&nextbit,outbuf);
  prefix = ic[0];
  /*-----------------loop over all stuff--------------*/
  for (i=1; i<imax; i++)
  {  now=ic[i];
    /*----------------search table----------------------*/
    index=-1; 
    jstart=prefix-firstcode; if (jstart<0) jstart=0;      
    for (j=jstart;j<itable;j++)
    {       
      if (now==addchar[j] && prefix==prefixcode[j])
      { index = j; prefix = lzwcode[j];
        break;   }
    }
    /*---------------add to table and dump code ---------*/
    if (index==-1)
    { p_writebits(prefix, nbits,fp,&size,&nextbit,outbuf);
      if ( nextcode == (1<<nbits)) 
      {   nbits++;
        if (nbits>maxbits)  /* full code list, reset table */
        { p_writebits(clearcode, nbits-1,fp,&size,&nextbit,outbuf);
          itable=0; nbits=jtot+1; nextcode=firstcode;
          prefix=now; i++; now=ic[i];
          p_writebits(prefix, nbits,fp,&size,&nextbit,outbuf);
        }
      }
      prefixcode[itable] = prefix;
      addchar[itable]=now;
      lzwcode[itable]=nextcode;
      nextcode++; 
      itable++; 
      prefix = now;
    }
  }
  /*---------------last character------------------*/
  p_writebits(prefix, nbits,fp,&size,&nextbit,outbuf);
  p_writebits(endcode, nbits,fp,&size,&nextbit,outbuf);
  p_writebits(0, -1,fp,&size,&nextbit,outbuf);
}

