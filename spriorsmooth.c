// Sphere prior shape model based smoothing of polygon meshes
// Copyright (C) Jussi Tohka, Institute of Signal Processing, Tampere University of
// Technology, 2004 - 2007.
// P.O. Box 553, FIN-33101, Finland
//  E-mail: jussi.tohka@tut.fi
// ----------------------------------------------------------------------------
// Permission to use, copy, modify, and distribute this software 
// for any purpose and without fee is hereby
// granted, provided that the above copyright notice appear in all
// copies.  The author and Tampere University of Technology make no representations
// about the suitability of this software for any purpose.  It is
// provided "as is" without express or implied warranty.
//  ----------------------------------------------------------------------------

// -----------------------------------------------------------------
// The method is described in 
// J. Tohka. 
// Surface Smoothing Based on a Sphere Shape Model. 
// In Proc. of 6th Nordic Signal Processing Symposium, NORSIG2004, 
// pp. 17 - 20, Espoo, Finland, June, 2004.
// ----------------------------------------------------------------




#include <stdio.h>
#include <stdlib.h>
#include "ssmooth.h"


int main(int argc, char** argv)
{
   trimesh mesh;
   simplexmesh smesh;
   int error_status,tmp;
   double lambda,multiplier;   

   if(argc < 3) {
     printf(" \n Usage: spriorsmooth infile outfile [lambda] [multiplier] \n ");
     printf(" the parameter lambda specifies the amount of smoothing (default 400) \n"); 
     printf(" the parameter multiplier specifies the desired accuracy of computation (default 0.1). \n"); 
     printf(" for a speed up increment the value. \n"); 
     return(1);
   }
   if(argc > 3) {   
     lambda = atof(argv[3]);
   }
   else {
     lambda = 400;
   }
   if(argc > 4) {
     multiplier = atof(argv[4]);
   }
   else {
     multiplier = 0.1;
   }
   
   printf("Reading the mesh to be smoothed \n");
   error_status  = ReadTriMesh(argv[1],&mesh);
   if(error_status != 0) {
     printf("Error %d while reading file \n",error_status); 
     return(2);
   }
   printf("Converting triangle mesh to simplex mesh \n"); 
   tmp = Tri2SimplexLinear(&mesh,&smesh);
   printf("Smoothing \n");
   tmp = SphereSmooth(&smesh, lambda, 1000, multiplier);
   printf("Converting simplex mesh to traingle mesh \n");
   tmp = Simplex2TriGeom(&smesh,&mesh);   
   printf("Writing triangle mesh \n");
   error_status = WriteTriMesh(argv[2],&mesh);
   if(error_status != 0) {
     printf("Error %d while writing file \n",error_status); 
     return(3);
   }   
  
   return(0);
}
