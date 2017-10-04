/* contains iexpand */
/* expand an index into its 4 components */
void iexpand(int k, int *idim, int *i)
{
  int j;
  i[0]=k%idim[0]; j=k/idim[0];
  i[1]=j%idim[1]; j=j/idim[1];
  i[2]=j%idim[2];
  i[3]=j/idim[2];
}
