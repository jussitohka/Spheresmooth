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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ssmooth.h"


// Reads a triangular mesh. MNI .obj file (ASCII)
// Not very much checking is done so crashes may follow if incorrect files are used
// See http://www.bic.mni.mcgill.ca/users/david/FAQ/polygons_format.txt for the file format
// Note that this function does not allow different colors for different triangles.

int ReadTriMesh(char* filename,trimesh* input_mesh)
{
  FILE *meshfile;
  char tmpc; 
  int i;
   
  // mesh
  meshfile = fopen(filename,"r");
  if(meshfile == NULL) return(1);
  if(fscanf(meshfile,"%c",&tmpc) != 1) {
     fclose(meshfile);
     return(2);
  } 
  for(i = 0;i < 5;i++) {
    if(fscanf(meshfile,"%lf",&(input_mesh->props[i])) != 1) {
      fclose(meshfile);
      return(3);
    }
  }
  if(fscanf(meshfile,"%d",&input_mesh->vertices) != 1){
    fclose(meshfile);
    return(4);
  }
  input_mesh->x = malloc(input_mesh->vertices * sizeof(points));
  for(i = 0;i < input_mesh->vertices;i++) {
    fscanf(meshfile,"%lf",&input_mesh->x[i][0]);
    fscanf(meshfile,"%lf",&input_mesh->x[i][1]);
    fscanf(meshfile,"%lf",&input_mesh->x[i][2]);
  }
  input_mesh->normals = malloc(input_mesh->vertices * sizeof(points));
  for(i = 0;i < input_mesh->vertices;i++) {
    fscanf(meshfile,"%lf",&input_mesh->normals[i][0]);
    fscanf(meshfile,"%lf",&input_mesh->normals[i][1]);
    fscanf(meshfile,"%lf",&input_mesh->normals[i][2]);
  }
  if(fscanf(meshfile,"%d",&input_mesh->triangles) != 1) { 
    fclose(meshfile);
    return(5);
  }
  for(i = 0;i < 5;i++) {
    if(fscanf(meshfile,"%d",&(input_mesh->props2[i])) != 1) {
      fclose(meshfile);
      return(6);
    }
  }
  input_mesh->triprops = malloc(input_mesh->triangles * sizeof(int));
  for(i = 0;i < input_mesh->triangles;i++) {
    fscanf(meshfile,"%d",&input_mesh->triprops[i]);
  }
  input_mesh->tri = malloc(input_mesh->triangles * sizeof(intpoints));    
  for(i = 0;i < input_mesh->triangles;i++) {
    fscanf(meshfile,"%d",&input_mesh->tri[i][0]);
    fscanf(meshfile,"%d",&input_mesh->tri[i][1]);
    fscanf(meshfile,"%d",&input_mesh->tri[i][2]);
  }
  //  for(i = 0;i < input_mesh->triangles;i++) {
  //  input_mesh->tri[i][0] =  input_mesh->tri[i][0] - 1;
  //  input_mesh->tri[i][1] =  input_mesh->tri[i][1] - 1;
  //  input_mesh->tri[i][2] =  input_mesh->tri[i][2] - 1;
  // }
    
  fclose(meshfile);
  return(0);
}   

// writes a triangular mesh in obj format. Normals are not updated.
// overwrites existing files.

int WriteTriMesh(char* filename,trimesh* output_mesh)
{
  FILE *meshfile;
  int i;
  
  meshfile = fopen(filename,"w");
  if(meshfile == NULL) return(1);
  
  fprintf(meshfile,"P ");
  for(i = 0;i < 5;i++) {
    fprintf(meshfile,"%lf ",output_mesh->props[i]);
  }
  fprintf(meshfile,"\n");
  fprintf(meshfile,"%d \n",output_mesh->vertices);
  for(i = 0;i < output_mesh->vertices;i++) {
    fprintf(meshfile,"%lf %lf %lf \n",output_mesh->x[i][0],output_mesh->x[i][1],output_mesh->x[i][2]);
  }
  for(i = 0;i < output_mesh->vertices;i++) {
    fprintf(meshfile,"%lf %lf %lf \n",output_mesh->normals[i][0],
            output_mesh->normals[i][1],output_mesh->normals[i][2]);
  }
  fprintf(meshfile,"%d\n",output_mesh->triangles);
  for(i = 0;i < 5;i++) {
    fprintf(meshfile,"%d ",output_mesh->props2[i]);
  }
  fprintf(meshfile,"\n");
   for(i = 0;i < output_mesh->triangles;i++) {
    fprintf(meshfile,"%d\n",output_mesh->triprops[i]);
  }
  
  for(i = 0;i < output_mesh->triangles;i++) {
    fprintf(meshfile,"%d %d %d \n",output_mesh->tri[i][0],output_mesh->tri[i][1],output_mesh->tri[i][2]);
  }
 
  fclose(meshfile);
  return(0);
}

int Tri2Simplex(trimesh* input_mesh,simplexmesh* out_mesh)
{
  //  simplexmesh* out_mesh; // simplex mesh to output
  int i,j;
  int n1,n2,n3;
  int triangles_found;
  int no_elem_tmp;
 
 
  // first compute number of vertices
  out_mesh->vertices = input_mesh->triangles;
  // allocate memory
 
  out_mesh->x = malloc(out_mesh->vertices * sizeof(points));
  out_mesh->nerel = malloc(out_mesh->vertices * sizeof(intpoints));
 
  // find the coordinates of each vertex in the simplex mesh
  // and find the neighbors of each vertex
  for(i = 0;i < input_mesh->triangles;i++) {
    // abbreaviations
    
    n1 = input_mesh->tri[i][0];
    n2 = input_mesh->tri[i][1];
    n3 = input_mesh->tri[i][2];
    
    out_mesh->x[i][0] = MEAN3(input_mesh->x[n1][0],input_mesh->x[n2][0],input_mesh->x[n3][0]);
    out_mesh->x[i][1] = MEAN3(input_mesh->x[n1][1],input_mesh->x[n2][1],input_mesh->x[n3][1]);
    out_mesh->x[i][2] = MEAN3(input_mesh->x[n1][2],input_mesh->x[n2][2],input_mesh->x[n3][2]);
    
    // find those triangles which share two of these nodes
    triangles_found = 0;
    j = 0;
   
    while(triangles_found < 3) {    
      if(COMMON_NODES(n1,n2,n3,input_mesh->tri[j][0],input_mesh->tri[j][1],input_mesh->tri[j][2]) == 2) { 
        out_mesh->nerel[i][triangles_found] = j;
        triangles_found++;        
      }
      j++;
    }
  }  
  
  // then construct faces matrix. for our purposes it does not have to be ordered. 
  out_mesh->no_faces = input_mesh->vertices;
  out_mesh->faces = malloc(input_mesh->vertices * sizeof(face));
 
  for(i = 0;i < out_mesh->no_faces;i++) {
    out_mesh->faces[i].no_elem = 0;
  }
  for(i = 0;i < input_mesh->triangles;i++) {
    for(j = 0;j < 3;j++) {
      no_elem_tmp = out_mesh->faces[input_mesh->tri[i][j]].no_elem;
      if(no_elem_tmp > (MAX_FACE_ELEMENTS - 1)) {
        printf(" ERROR: MAX_FACE_ELEMENTS too small \n \n");
        return(1); 
      } 
      out_mesh->faces[input_mesh->tri[i][j]].elem[no_elem_tmp] = i;
      out_mesh->faces[input_mesh->tri[i][j]].no_elem++; 
    }  
  }
 
  return(0);
}

// Linear time version of the previous function. 
// Requires somewhat more storage than the previous function and is somewhat more complicated

int Tri2SimplexLinear(trimesh* input_mesh,simplexmesh* out_mesh)
{
 
  int i,j,k,l;
  int t1,t2;
  int already_known_neighbour;
  int n1,n2,n3;
  int no_elem_tmp;
  // face* vertex_nerel;
  int* nerel_no;

  // first compute number of vertices
  out_mesh->vertices = input_mesh->triangles;
  // allocate memory
  out_mesh->x = malloc(out_mesh->vertices * sizeof(points));
  out_mesh->nerel = malloc(out_mesh->vertices * sizeof(intpoints));

 
  // STEP 1: construct faces matrix. for our purposes it does not have to be ordered. 
   printf("Step1 \n");
  out_mesh->no_faces = input_mesh->vertices;
  out_mesh->faces = malloc(input_mesh->vertices * sizeof(face));
  for(i = 0;i < out_mesh->no_faces;i++) {
    out_mesh->faces[i].no_elem = 0;
  }
 
  for(i = 0;i < input_mesh->triangles;i++) {
    for(j = 0;j < 3;j++) {
      no_elem_tmp = out_mesh->faces[input_mesh->tri[i][j]].no_elem;
      if(no_elem_tmp > (MAX_FACE_ELEMENTS - 1)) {
        printf(" ERROR: MAX_FACE_ELEMENTS too small \n \n");
        return(1); 
      } 
      out_mesh->faces[input_mesh->tri[i][j]].elem[no_elem_tmp] = i;
      out_mesh->faces[input_mesh->tri[i][j]].no_elem++; 
    }  
  }
  
  // STEP 2: then find the neighborhood relations.
  // For this we utilize the information about the faces matrix
 
  nerel_no = malloc(out_mesh->vertices * sizeof(int)); // the number of neighbours already 
                                                       //collected is stored in this array
  for(i = 0;i < out_mesh->vertices;i++) {
    nerel_no[i] = 0;
  }    

  for(i = 0;i < out_mesh->no_faces; i++) {
    for(j = 0;j < out_mesh->faces[i].no_elem;j++) {
      for(k = (j + 1);k < out_mesh->faces[i].no_elem;k++) {
        t1 = out_mesh->faces[i].elem[j];  // abbreviation
        t2 = out_mesh->faces[i].elem[k];  // abbreviation
        // if the next if statement is TRUE, t1 and t2 are neighbours
        if((COMMON_NODES(input_mesh->tri[t1][0],input_mesh->tri[t1][1],input_mesh->tri[t1][2],
                         input_mesh->tri[t2][0],input_mesh->tri[t2][1],input_mesh->tri[t2][2]) == 2)) {
          // then we figure out is this neighbourhood relation already known or not.
          // Here, there are possibilities to speed up the algorithm.  
          already_known_neighbour = 0;
          for(l = 0;l < nerel_no[t1];l++) {
            if(out_mesh->nerel[t1][l] == t2) already_known_neighbour = 1;
          }
          if(!already_known_neighbour) { 
            out_mesh->nerel[t1][nerel_no[t1]] = t2;
            out_mesh->nerel[t2][nerel_no[t2]] = t1;
            nerel_no[t1]++;
            nerel_no[t2]++;
          }
        }
      }
    }
  } 
  printf("Step 2 \n");

  // STEP3: Finally, find the coordinates of each vertex in the simplex mesh
  for(i = 0;i < input_mesh->triangles;i++) {
    // abbreaviations
    n1 = input_mesh->tri[i][0];
    n2 = input_mesh->tri[i][1];
    n3 = input_mesh->tri[i][2];
    out_mesh->x[i][0] = MEAN3(input_mesh->x[n1][0],input_mesh->x[n2][0],input_mesh->x[n3][0]);
    out_mesh->x[i][1] = MEAN3(input_mesh->x[n1][1],input_mesh->x[n2][1],input_mesh->x[n3][1]);
    out_mesh->x[i][2] = MEAN3(input_mesh->x[n1][2],input_mesh->x[n2][2],input_mesh->x[n3][2]);
  }  
  return(0);
}


// this function is not required for surface smoothing
void RemoveDoubleNeighbours(face* pface) 
{
  int i,j;
  face tmp_face;
  int is_in_face;  

  // we utilize the facts that the first two elements are not doubles 
  // and the last two are surely doubles to speed up the computation
  tmp_face.no_elem = 2;
  tmp_face.elem[0] = pface->elem[0];
  tmp_face.elem[1] = pface->elem[1];
  for(i = 2; i < (pface->no_elem - 2);i++) {
    is_in_face = 0;
    for(j = 0;j < i;j++) {
      if(tmp_face.elem[j] == pface->elem[i]) is_in_face = 1;
    }
    if(!is_in_face) {
      tmp_face.elem[tmp_face.no_elem] = pface->elem[i];
      tmp_face.no_elem++; 
    }
  }
  pface->no_elem = tmp_face.no_elem;
  for(i = 2;i < pface->no_elem;i++) {
    pface->elem[i] = tmp_face.elem[i];
  }
}


// this function is not required for surface smoothing

trimesh* Simplex2Tri(simplexmesh* input_mesh) 
{
}

int CommonNodes(int n1, int n2, int n3, int m1, int m2, int m3)
{
  int sum = 0;
  if( (n1 == m1) || (n1 == m2) || (n1 == m3)) sum++;
  if( (n2 == m1) || (n2 == m2) || (n2 == m3)) sum++;
  if( (n3 == m1) || (n3 == m2) || (n3 == m3)) sum++;
  return(sum);
}

// The main thing: sphere smoothing algorithm for simplex meshes

int SphereSmooth(simplexmesh* input_mesh,double lambda,int max_iter, double multiplier) 
{
  int iter_counter = 0;
  double alpha,alpha2,eta,theta,theta_tmp,threshold;
  points* x_old;
  points* x_new;  
  points* ip;
  points* ipsqr;
  points g;
  int i,j;
  int n1,n2,n3; 

  alpha = 1/(3*cos(2*atan((2 * sqrt(PI * sqrt(3)))/(3 * sqrt(input_mesh->vertices)))));
  alpha2 = (3*alpha - 1)*(3*alpha - 1)/input_mesh->vertices; 
  eta = 0.98*(2/(1 + lambda*(1 + 3*alpha)*(1 + 3*alpha))); // this guarantees convergence. 
                                                           // Proving this turned out be quite easy matrix
                                                           // algebra exercise
  threshold = ComputeThreshold(input_mesh, multiplier);
  theta = threshold + 1;  
  
  x_old = malloc(input_mesh->vertices * sizeof(points));
  x_new = malloc(input_mesh->vertices * sizeof(points));
  ip =  malloc(input_mesh->vertices * sizeof(points));
  ipsqr =  malloc(input_mesh->vertices * sizeof(points));

  for(i = 0;i < input_mesh->vertices;i++) {
    for(j = 0;j < 3;j++) {
      x_new[i][j] = input_mesh->x[i][j];
    } 
  }

  while((iter_counter < max_iter) & (theta > threshold)) {
    for(i = 0;i < input_mesh->vertices;i++) {
      for(j = 0;j < 3;j++) {
        x_old[i][j] = x_new[i][j];
      } 
    }
    for(i = 0;i < input_mesh->vertices;i++) {
      // abbreavations
      n1 = input_mesh->nerel[i][0];
      n2 = input_mesh->nerel[i][1];
      n3 = input_mesh->nerel[i][2];
      for(j = 0;j < 3;j++) {
        ip[i][j] = x_old[i][j] - alpha*(x_old[n1][j] + x_old[n2][j] + x_old[n3][j]);
      }  
    }  
    for(i = 0;i < input_mesh->vertices;i++) {
      // abbreavations
      n1 = input_mesh->nerel[i][0];
      n2 = input_mesh->nerel[i][1];
      n3 = input_mesh->nerel[i][2];
      for(j = 0;j < 3;j++) {
        ipsqr[i][j] = ip[i][j] - alpha*(ip[n1][j] + ip[n2][j] + ip[n3][j]);
      }  
    }
    g[0] = 0;
    g[1] = 0;
    g[2] = 0;  
    for(i = 0;i < input_mesh->vertices;i++) {     
      for(j = 0;j < 3;j++) {
        g[j] = g[j] + alpha2*x_old[i][j];
      }
    }
    for(i = 0;i < input_mesh->vertices;i++) {     
      for(j = 0;j < 3;j++) {
        x_new[i][j] = x_old[i][j] - eta*(x_old[i][j] + lambda * ipsqr[i][j] + g[j]) + eta*(input_mesh->x[i][j]);
      } 
    }
    theta = 0;
    for(i = 0;i < input_mesh->vertices;i++) {
      for(j = 0;j < 3;j++) {
         theta_tmp = fabs(x_new[i][j] - x_old[i][j]);
         if(theta_tmp > theta) theta = theta_tmp;
      }
    }
    iter_counter++;
    printf("*");  
    
  }
  printf("\n");

  for(i = 0;i < input_mesh->vertices;i++) {
    for(j = 0;j < 3;j++) {
      input_mesh->x[i][j] = x_new[i][j]; 
    }
  } 
  free(x_new);
  free(x_old);
  free(ip);
  free(ipsqr); 
  return(0);
}

double ComputeThreshold(simplexmesh* input_mesh, double multiplier)
{

  double threshold = 0;
  double threshold2 = 0;
  int i,j;
  points g;
 
  g[0] = 0;
  g[1] = 0;
  g[2] = 0;
  for(i = 0;i < input_mesh->vertices;i++) {     
    for(j = 0;j < 3;j++) {
        g[j] = g[j] + (input_mesh->x[i][j])/((double) input_mesh->vertices);
    }
  }
  for(i = 0;i < input_mesh->vertices;i++) {
    threshold2 = 0;
    for(j = 0;j < 3;j++) {
       threshold2 = threshold2 + (input_mesh->x[i][j] - g[j])*(input_mesh->x[i][j] - g[j]);
    }
    threshold = threshold + sqrt(threshold2)/(((double) input_mesh->vertices)*((double) input_mesh->vertices));
  }
  threshold = multiplier*threshold;
  return(threshold); 
}

int Simplex2TriGeom(simplexmesh* input_smesh,trimesh* input_trimesh)
{
  int i,j,k,nf;
  double coor;
   
  for(i = 0;i < input_trimesh->vertices;i++) {
    // 
    nf = input_smesh->faces[i].no_elem;
    for(j = 0;j < 3;j++) { 
      coor = 0;
      for(k = 0;k <nf;k++) {
        coor = coor + input_smesh->x[input_smesh->faces[i].elem[k]][j];
      }   
      coor = coor / ((double) nf);
      input_trimesh->x[i][j] = coor;
    }
  }
  return(0);
}
