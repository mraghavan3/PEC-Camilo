#include "velocityMeshUtilities.h"
#define TYPE_TET 0
#define TYPE_HEX 1
#define TYPE_TRI 2

struct Cell{

	vector<int> elementID;

};

int maxNumber (int a, int b){
	return a>b? a:b;
}

int minNumber (int a, int b){
	return a<b? a:b;   
}

bool testPointandTri(double x, double y, double z, 
					double x1, double y1, double z1, 
					double x2, double y2, double z2, 
					double x3, double y3, double z3,bool printDebug
					){

						/*cout << "x " << x << " y " << y 
							<< " x1 " << x1 << " y1 " << y1 
							<< " x2 " << x2 << " y2 " << y2 
							<< " x3 " << x3 << " y3 " << y3<<endl; */
/*double x12 = x1-x2;*/ double x13 = x1-x3; 
double x21 = x2-x1 ; //double x23 = x2-x3; 
/*double x31 = x3-x1;*/ double x32 = x3-x2; 

double y12 = y1-y2; //double y13 = y1-y3; 
/*double y21 = y2-y1 ;*/ double y23 = y2-y3; 
double y31 = y3-y1; //double y32 = y3-y2; 
/*
cout << "x13 " << x13 << " x21 " << x21 
							<< " x32 " << x32 << " y12 " << y12 
							<< " y23 " << y23 << " y31 " << y2 
							<<endl;*/

// a tri is in 2D

/*
double z12 = z1-z2; double z13 = z1-z3; 
double z21 = z2-z1 ; double z23 = z2-z3; 
double z31 = z3-z1; double z32 = z3-z2; 
*/
double A = (x2*y3-x3*y2)+(x3*y1-x1*y3)+(x1*y2-x2*y1);


double A23 = x2*y3-x3*y2; 
double A31 = x3*y1-x1*y3;
double A12 = x1*y2-x2*y1;
//cout << "A23 " << A23 << " A31 " << A31 <<" A12 " << A12 << " A " << A <<endl;
double Xi1 = (A23+ y23*x + x32*y  ) /A  ; 
double Xi2 = (A31 +y31*x + x13*y  ) /A  ;
double Xi3 = (A12 + y12*x + x21*y  ) /A  ; 

if(printDebug){
	cout << "x " << x << " y " << y 
							<< " x1 " << x1 << " y1 " << y1 
							<< " x2 " << x2 << " y2 " << y2 
							<< " x3 " << x3 << " y3 " << y3<<endl; 
	cout << "x13 " << x13 << " x21 " << x21 
							<< " x32 " << x32 << " y12 " << y12 
							<< " y23 " << y23 << " y31 " << y2 
							<<endl;
cout << "A23 " << A23 << " A31 " << A31 <<" A12 " << A12 << " A " << A <<endl;
cout<<" Xi1 " <<Xi1<<" Xi2 " << Xi2 << " Xi3 "<<Xi3<<endl;
}

//cout<< " solution " << A12 + y12*x + x21*y <<" ";
//cout<<" Xi1 " <<Xi1<<" Xi2 " << Xi2 << " Xi3 "<<Xi3<<endl;

if(abs(Xi1)<0.0001) {Xi2-=abs(Xi1);  Xi1=0; } 
if(abs(Xi2)<0.0001) {Xi1-=abs(Xi2);Xi2=0; } 
if(abs(Xi3)<0.0001) {Xi1-=abs(Xi3);Xi3=0; } 
return (Xi1>=0 && Xi2>=0 && Xi3>=0 && abs((Xi1+Xi2+Xi3 )-1) <0.0001 );

}


void checkHingesElements (const vector<double> &velMeshPoints, const vector<int>&velMeshConnectivity,const vector<int> meshNeighbors,vector<fiber> &fibers, int numElements, int numnodes,int meshType, int nodesPerElement ){

	//cout<<" fiber 0 hinge 0 inside function" <<  fibers[0].hinges[0].position(0)<<endl;

	for (int fiberID = 0; fiberID < fibers.size(); fiberID++)
	{
		for (int hingeID = 0; hingeID < fibers[fiberID].numberOfHinges; hingeID++)
		{
			

			double hingeX = fibers[fiberID].hinges[hingeID].position[0];
			double hingeY = fibers[fiberID].hinges[hingeID].position(1);
			double hingeZ = fibers[fiberID].hinges[hingeID].position(2);
			int elementID = fibers[fiberID].hinges[hingeID].elementID;

			if( !fibers[fiberID].hinges[hingeID].isStationary ){
		 // cout<< "fiber " << fiberID << " Hinge " << hingeID<<" elementID" <<elementID << " x "<< hingeX << " y" << hingeY << endl;
			//for tri elements get rid of z corods

			if ( testPointandTri(hingeX,hingeY,0 ,
				velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+0] + numnodes * 0     ], //first node(0) position x(0) 
				velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+0] + numnodes * 1     ], //first node(0) position y(1) 
				0, //first node(0) position z(2) 
				velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+1] + numnodes * 0     ], //second node(1) position x(0) 
				velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+1] + numnodes * 1     ], //second node(1) position y(1) 
				0, //second node(1) position z(2) 
				velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+2] + numnodes * 0     ], //third node(2) position x(0) 
				velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+2] + numnodes * 1     ], //third node(2) position y(1) 
				0 //third node(2) position z(2) 
			,false)){
				// the hinge is in the triangle

			}else{


				// look in the neighboring triangles

				for (int neighbor = 0; neighbor < 3; neighbor++)
				{
					int neighborID = meshNeighbors[elementID*3+neighbor ];

					if (neighborID!=-1)
					{
						
						//cout<<"size "<<  velMeshConnectivity.size()<< " index : " <<  neighborID*6+0 + numnodes * 0;
						
						//for tri elements get rid of z corods
						if ( testPointandTri(hingeX,hingeY,0,
							velMeshPoints[  velMeshConnectivity[neighborID*nodesPerElement+0] + numnodes * 0     ], //first node(0) position x(0) 
							velMeshPoints[  velMeshConnectivity[neighborID*nodesPerElement+0] + numnodes * 1     ], //first node(0) position y(1) 
							0, //first node(0) position z(2) 
							velMeshPoints[  velMeshConnectivity[neighborID*nodesPerElement+1] + numnodes * 0     ], //second node(1) position x(0) 
							velMeshPoints[  velMeshConnectivity[neighborID*nodesPerElement+1] + numnodes * 1     ], //second node(1) position y(1) 
							0, //second node(1) position z(2) 
							velMeshPoints[  velMeshConnectivity[neighborID*nodesPerElement+2] + numnodes * 0     ], //third node(1) position x(0) 
							velMeshPoints[  velMeshConnectivity[neighborID*nodesPerElement+2] + numnodes * 1     ], //third node(1) position y(1) 
							0 //third node(2) position z(2) 
						,false)){ 
							//element updated
							fibers[fiberID].hinges[hingeID].elementID=neighborID;
							//cout<<" element found "<<endl; 
 						}
						else{ 
							//element was not one of the neighbors, have to check again. 
							//cout<<" fiber "<< fiberID << " hinge " <<  hingeID << " does not have an associated element to it, the last element was"<< elementID  <<endl;
						}


					}
				}
			}
			}
		}//for
	}//for
}




bool isPointInElement(double pointX, double pointY, double pointZ, vector<double> meshPoints, vector<int> meshConnectivity, int elementID, int elementType    ){

	if (elementType == TYPE_TRI)
	{



	}


}


void findElementMidpoint(vector<double> &nodesInfo, vector<int> &meshConnectivity,int numElements, int numNodes, int elemntID, double &midX, double& midy, double & midz, int elementType, int nodesPerElement){

	double minx,miny,minz,maxx,maxy,maxz;

	if (elementType == TYPE_HEX)
	{
		// if it is a type hex it must have 8 nodes. 
		// get the minx, maxx, miny, maxy, minz, maxz
		// the mid point x will be (maxx + minx / 2),
		// same for all other coordinates
		
		int coordStride = numNodes*numElements;
		
		// node 0 of current element
		minx = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*0]];
		maxx = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*0]];
		miny = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*1]];
		maxy = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*1]];
		minz = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*2]];
		maxz = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*2]];

		for (int nodeID = 1; nodeID < numNodes; nodeID++)
		{
			if( minx > nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*0 + nodeID  ]]   ) 
			{ minx = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*0 + nodeID  ]]; }

			if( maxx < nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*0 + nodeID  ]]   ) 
			{ maxx = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*0 + nodeID  ]]; }

			if( miny > nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*1 + nodeID  ]]   ) 
			{ miny = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*1 + nodeID  ]]; }

			if( maxy < nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*1 + nodeID  ]]   ) 
			{ maxy = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*1 + nodeID  ]]; }

			if( minz > nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*2 + nodeID  ]]   ) 
			{ minz = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*2 + nodeID  ]]; }

			if( maxz < nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*2 + nodeID  ]]   ) 
			{ maxz = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*2 + nodeID  ]]; }

			midX = (minx+maxx)/ 2;
			midy = (miny+maxy)/ 2;
			midz = (minz+maxz)/ 2;
		}
	}


	if (elementType == TYPE_TRI)
	{
		// if it is a type tri it must have 3 nodes. 
		// get the minx, maxx, miny, maxy, minz, maxz
		// the mid point x will be (maxx + minx / 2),
		// same for all other coordinates
		
		int coordStride = numNodes; // todo change this
		
		
		// node 0 of current element
		minx = nodesInfo[meshConnectivity[elemntID*nodesPerElement] ];
		maxx = nodesInfo[meshConnectivity[elemntID*nodesPerElement] ];
		miny = nodesInfo[meshConnectivity[elemntID*nodesPerElement] + numNodes ];
		maxy = nodesInfo[meshConnectivity[elemntID*nodesPerElement] + numNodes ];
		
		for (int nodeID = 1; nodeID < nodesPerElement; nodeID++)
		{
			if( minx > nodesInfo[meshConnectivity[elemntID*nodesPerElement +  nodeID  ]]   ) 
			{ minx = nodesInfo[meshConnectivity[elemntID*nodesPerElement +  nodeID  ]]; }

			if( maxx < nodesInfo[meshConnectivity[elemntID*nodesPerElement + nodeID  ]]   ) 
			{ maxx = nodesInfo[meshConnectivity[elemntID*nodesPerElement +  nodeID  ]]; }

			if( miny > nodesInfo[meshConnectivity[elemntID*nodesPerElement  + nodeID  ] + coordStride*1 ]   ) 
			{ miny = nodesInfo[meshConnectivity[elemntID*nodesPerElement  + nodeID  ]+ coordStride*1]; }

			if( maxy < nodesInfo[meshConnectivity[elemntID*nodesPerElement  + nodeID  ]+ coordStride*1]   ) 
			{ maxy = nodesInfo[meshConnectivity[elemntID*nodesPerElement + nodeID  ]+ coordStride*1 ]; }

			
			
		}
		midX = (minx+maxx)/ 2;
			midy = (miny+maxy)/ 2;
			
		//cout<< "mid x " << midX << " mid y " << midy<< " min x " << minx<<" min y " << miny<<  "maxx"<< maxx<< "maxy" << maxy <<endl  ;
	}
}

void findElementsize(vector<double> &nodesInfo, vector<int> &meshConnectivity,int numElements, int numNodes, int elemntID, double &sizex, double& sizey, double & sizez, int elementType, int nodesPerElement){

	double minx,miny,minz,maxx,maxy,maxz;

	if (elementType == TYPE_HEX)
	{
		// if it is a type hex it must have 8 nodes. 
		// get the minx, maxx, miny, maxy, minz, maxz
		// the mid point x will be (maxx + minx / 2),
		// same for all other coordinates
		
		int coordStride = numNodes*numElements;
		
		// node 0 of current element
		minx = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*0]];
		maxx = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*0]];
		miny = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*1]];
		maxy = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*1]];
		minz = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*2]];
		maxz = nodesInfo[meshConnectivity[elemntID*numNodes+coordStride*2]];

		for (int nodeID = 1; nodeID < numNodes; nodeID++)
		{
			if( minx > nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*0 + nodeID  ]]   ) 
			{ minx = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*0 + nodeID  ]]; }

			if( maxx < nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*0 + nodeID  ]]   ) 
			{ maxx = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*0 + nodeID  ]]; }

			if( miny > nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*1 + nodeID  ]]   ) 
			{ miny = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*1 + nodeID  ]]; }

			if( maxy < nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*1 + nodeID  ]]   ) 
			{ maxy = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*1 + nodeID  ]]; }

			if( minz > nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*2 + nodeID  ]]   ) 
			{ minz = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*2 + nodeID  ]]; }

			if( maxz < nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*2 + nodeID  ]]   ) 
			{ maxz = nodesInfo[meshConnectivity[elemntID*numNodes + coordStride*2 + nodeID  ]]; }

			
		}

		sizex = (minx-maxx);
			sizey = (miny-maxy);
			sizez = (minz-maxz);
			

	}


	if (elementType == TYPE_TRI)
	{
		// if it is a type tri it must have 3 nodes. 
		// get the minx, maxx, miny, maxy, minz, maxz
		// the mid point x will be (maxx + minx / 2),
		// same for all other coordinates
		
		int coordStride = numNodes;
		
		// node 0 of current element
		minx = nodesInfo[meshConnectivity[elemntID*nodesPerElement] ];
		maxx = nodesInfo[meshConnectivity[elemntID*nodesPerElement] ];
		miny = nodesInfo[meshConnectivity[elemntID*nodesPerElement] + numNodes ];
		maxy = nodesInfo[meshConnectivity[elemntID*nodesPerElement] + numNodes ];
		
		for (int nodeID = 1; nodeID < nodesPerElement; nodeID++)
		{
			if( minx > nodesInfo[meshConnectivity[elemntID*nodesPerElement +  nodeID  ]]   ) 
			{ minx = nodesInfo[meshConnectivity[elemntID*nodesPerElement +  nodeID  ]]; }

			if( maxx <nodesInfo[meshConnectivity[elemntID*nodesPerElement + nodeID  ]]   ) 
			{ maxx = nodesInfo[meshConnectivity[elemntID*nodesPerElement +  nodeID  ]]; }

			if( miny > nodesInfo[meshConnectivity[elemntID*nodesPerElement  + nodeID  ] + coordStride*1 ]   ) 
			{ miny = nodesInfo[meshConnectivity[elemntID*nodesPerElement  + nodeID  ]+ coordStride*1]; }

			if( maxy < nodesInfo[meshConnectivity[elemntID*nodesPerElement  + nodeID  ]+ coordStride*1]   ) 
			{ maxy = nodesInfo[meshConnectivity[elemntID*nodesPerElement + nodeID  ]+ coordStride*1 ]; }
		}

		   sizex = abs(minx-maxx);
			sizey = abs(miny-maxy);
			sizez = abs(minz-maxz);

			//cout<< " sizex " << sizex << " sizey " << sizey <<" sizez " << sizez<<endl;
			//cout<< " minx " << minx << " maxx " << maxx <<" miny " << miny<<" maxy " << maxy<<endl;
	}
}



void sort_heap_external (int n, int &indx , int & i, int &j, int isgn){

  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( indx < 0 )
  {
    if ( indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      indx = 2;
      i = i_save;
      j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        indx = 1;
      }
      i = i_save;
      j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( indx == 1 ) 
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      indx = -2;
      i = i_save;
      j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    indx = 0;
    i = i_save;
    j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    indx = 1;
    i = i_save;
    j = j_save;
  }

  return;
}


void i4col_swap ( int m, int n, vector<int> &matrixA, int icol1, int icol2 ){

# define OFFSET 1

  int i;
  int t;
//
//  Check.
//
  if ( icol1 - OFFSET < 0 || n-1 < icol1 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL1 is out of range.\n";
    exit ( 1 );
  }

  if ( icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL2 is out of range.\n";
    exit ( 1 );
  }

  if ( icol1 == icol2 )
  {
    return;
  }
  for ( i = 0; i < m; i++ )
  {
    t                     = matrixA[i+(icol1-OFFSET)*m];
    matrixA[i+(icol1-OFFSET)*m] = matrixA[i+(icol2-OFFSET)*m];
    matrixA[i+(icol2-OFFSET)*m] = t;
  }

  return;
# undef OFFSET

}

int i4col_compare ( int m, int n, vector<int> &matrixA, int i, int j ){

 int k;
//
//  Check.
//
  if ( i < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index I = " << i << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < i )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index I = " << i << ".\n";
    exit ( 1 );
  }

  if ( j < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index J = " << j << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < j )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index J = " << j << ".\n";
    exit ( 1 );
  }

  if ( i == j )
  {
    return 0;
  }

  k = 1;

  while ( k <= m )
  {
	  if ( matrixA[k-1+(i-1)*m] < matrixA[k-1+(j-1)*m] )
    {
      return (-1);
    }
    else if ( matrixA[k-1+(j-1)*m] < matrixA[k-1+(i-1)*m] )
    {
      return 1;
    }
    k = k + 1;
  }

  return 0;
}


void sortColumnsofMatrix(int rows, int columns, vector<int> &arrayA){

	//taken from http://people.sc.fsu.edu/~jburkardt/cpp_src/tet_mesh_tet_neighbors/tet_mesh_tet_neighbors.html
	int i,indx,isgn, j;


	//initialize
	i =0;
	indx = 0;
	isgn =0;
	j =0;

	// use external heap sorter

	for( ; ; )
	{
		sort_heap_external(columns,indx,i,j,isgn );

		if ( 0 < indx )
		{
			i4col_swap(rows,columns,arrayA,i,j);
		}

		else if ( indx<0 )
		{
			isgn = i4col_compare ( rows, columns, arrayA, i, j );
		}
		else if ( indx==0 )
		{
			break;
		}

	}
	return;
}


int find_tri_neighbors(int numTris, int nodesPerTri, vector<int> & tri_points, vector<int> & triNeighbors){

	vector<int> faces;
	int numColumns = 4;
	int numFacesPerTri = 3;

	int triStride = numColumns*numFacesPerTri;
	int faceStride = numColumns;

	faces.resize(numColumns*numFacesPerTri*numTris); 
    // the edges of a triangle with index T and nodes I,J,K, will be stored as follows:
	// (I,J,0,T, j,K,1,T, K,i,2,T)
	// Thus there are 4 "columns" and 3 faces per tri
	// Also the indices I,J,K  will be sorted to make matching easier afterwards

	int a,b,i,j,k,faceID;


	for (int tri = 0; tri < numTris; tri++)
	{
		// pull node indices for current tet

		i = tri_points[tri*nodesPerTri+0];
		j = tri_points[tri*nodesPerTri+1];
		k = tri_points[tri*nodesPerTri+2];
		
		//first edge (I,J), face 0
		faceID =0;
		
		faces[tri*triStride + faceID*faceStride + 0] = minNumber(i,j);
		faces[tri*triStride + faceID*faceStride + 1] = maxNumber(i,j);
		faces[tri*triStride + faceID*faceStride + 2] = faceID;
		faces[tri*triStride + faceID*faceStride + 3] = tri;

		//second face (I,K), face 0
		faceID =1;
		
		faces[tri*triStride + faceID*faceStride + 0] = minNumber(j,k);
		faces[tri*triStride + faceID*faceStride + 1] = maxNumber(j,k);
		faces[tri*triStride + faceID*faceStride + 2] = faceID;
		faces[tri*triStride + faceID*faceStride + 3] = tri;

		//third face (J,K), face 0
		faceID =2;
		
		faces[tri*triStride + faceID*faceStride + 0] = minNumber(i,k);
		faces[tri*triStride + faceID*faceStride + 1] = maxNumber(i,k);
		faces[tri*triStride + faceID*faceStride + 2] = faceID;
		faces[tri*triStride + faceID*faceStride + 3] = tri;

	}

	
	/*
	cout<<" face size " << faces.size()<<endl;

	

	for (int pp = 0; pp < 12; pp++)
	{

	cout<<" tri " << floor(pp/3)  << " FaceID "<<pp%3<<endl;

		for (int qq = 0; qq < 2; qq++)
		{
		
			cout<<"index in array" << (floor(pp/3) *triStride + pp%3*faceStride + qq ) ;
			cout<< " node index " <<faces[floor(pp/3) *triStride + pp%3*faceStride + qq ]<<"  ,";
		
		}
	cout<<endl;
	}

	cout<<endl;

	*/


	sortColumnsofMatrix(faceStride,numFacesPerTri*numTris,faces);





	int face,face1,face2,tetra1,tetra2;

	for ( j = 0; j < numTris; j++ )
  {
    for ( i = 0; i < numFacesPerTri; i++ )
    {
		triNeighbors[i+j*numFacesPerTri] = -1;
    }
  }

  face = 0;

  for ( ; ; )
  {
    if ( numFacesPerTri * numTris - 1 <= face )
    {
      break;
    }

    if ( faces[0+face*numColumns] == faces[0+(face+1)*numColumns] &&
         faces[1+face*numColumns] == faces[1+(face+1)*numColumns]  )
    {
      face1 = faces[2+face*numColumns];
      tetra1 = faces[3+face*numColumns];
      face2 = faces[2+(face+1)*numColumns];
      tetra2 = faces[3+(face+1)*numColumns];

	  if(tetra1<numTris){ triNeighbors[face1+tetra1*3] = tetra2 ;}
	  if(tetra1<numTris){ triNeighbors[face2+tetra2*3] = tetra1 ;}
      face = face + 2;
    }
    else
    {
      face = face + 1;
    }
  }


	return 0;

}


void find_Initial_tetID(vector<double>&velMeshPoints,vector<int>&velMeshConnectivity,vector<fiber> &fibers, int numElements, int numnodes,int meshType,int nodesPerElement ){

	// preprocessing
	// find the enclosure of the fibers, that way we will reduce the number of elements we will have to store
	double minx,miny,minz,maxx, maxy,maxz;

	minx = fibers[0].hinges[0].position[0];
	miny = fibers[0].hinges[0].position[1];
	minz = fibers[0].hinges[0].position[2];

	maxx = fibers[0].hinges[0].position[0];
	maxy = fibers[0].hinges[0].position[1];
	maxz = fibers[0].hinges[0].position[2];

	for (int ii = 0; ii < fibers.size(); ii++)
	{
		for (int jj = 0; jj < fibers[ii].hinges.size(); jj++)
		{
			//cout<< " hinge x " << fibers[ii].hinges[jj].position[0] << " y " <<fibers[ii].hinges[jj].position[1] <<endl;
			if (fibers[ii].hinges[jj].position[0] > maxx){ maxx = fibers[ii].hinges[jj].position[0];	}
			if (fibers[ii].hinges[jj].position[0] < minx){ minx = fibers[ii].hinges[jj].position[0];	}

			if (fibers[ii].hinges[jj].position[1] > maxy){ maxy = fibers[ii].hinges[jj].position[1];	}
			if (fibers[ii].hinges[jj].position[1] < miny){ miny = fibers[ii].hinges[jj].position[1];	}
			
			if (fibers[ii].hinges[jj].position[2] > maxz){ maxz = fibers[ii].hinges[jj].position[2];	}
			if (fibers[ii].hinges[jj].position[2] < minz){ minz = fibers[ii].hinges[jj].position[2];	}

		}

	}
	vector<double> boxSize;
	boxSize.resize(3);

	
	// The grid will be as big as the dimensions of the mins and maxs,
	//ideally the grid should have the sime size as the mesh, for that we should compute the 
	//average element size in the 3 dimensions. That size is rc
	
	double sizex,sizey,sizez,maxsize ;
	maxsize =0;
	
	for (int element = 0; element < velMeshConnectivity.size()/nodesPerElement; element++) // 6 for quadratic triangles
	{
		findElementsize(velMeshPoints,velMeshConnectivity, numElements, numnodes, element,sizex,sizey,sizez,meshType,nodesPerElement);
		
		if (sizex>maxsize){	maxsize =sizex; }
		if (sizey>maxsize){	maxsize =sizey; }
		if (sizez>maxsize){	maxsize =sizez; }

	}


	double rc= maxsize;

	boxSize[0] =( maxx-minx)+rc;
	boxSize[1] =(maxy-miny)+rc;
	boxSize[2] = maxz-minz+rc;
	minx -= rc/2;
	miny -= rc/2;
	minz -= rc/2;

	//cout<<"box size " << boxSize[0] << " " <<boxSize[1]<<endl;
	
	vector<int> cellsNumber;
	cellsNumber.resize(3);

	cellsNumber[0] = (int) floor(boxSize[0] /rc);
	cellsNumber[1] = (int) floor(boxSize[1] /rc);
	cellsNumber[2] = (int) floor(boxSize[2] /rc);

	if (cellsNumber[0] == 0){cellsNumber[0] =1	;}
	if (cellsNumber[1] == 0){cellsNumber[1] =1	;}
	if (cellsNumber[2] == 0){cellsNumber[2] =1	;}
	int numCells = cellsNumber[0]* cellsNumber[1] * cellsNumber[2];

	//cout<< " cells number " << cellsNumber[0] << " " << cellsNumber[1]<< " " <<cellsNumber[2]<<endl;

	//cout<<" max size " << maxsize<<endl;

	vector<double> cellLength;
	cellLength.resize(3);
	cellLength[0] = boxSize[0]/cellsNumber[0];
	cellLength[1] = boxSize[1]/cellsNumber[1];
	cellLength[2] = boxSize[2]/cellsNumber[2];
	
	//cout<< " cells length " << cellLength[0] << " " << cellLength[1]<< " " <<cellLength[2]<<endl;
	
	vector<vector<int>> CellsVector; 
	
	

	for (int i = 0; i < numCells; i++)
	{
		vector<int> vec;
		CellsVector.push_back(vec);
	}
		

	// go through each element and locate each element center, if it is outside the grid do not count it

	double midx,midy,midz;
	int cellIDx, cellIDy, cellIDz,CellID;
	int elementsFound =0;
	for (int element = 0; element < velMeshConnectivity.size()/nodesPerElement; element++)
	{
		findElementMidpoint(velMeshPoints,velMeshConnectivity, numElements, numnodes, element,midx,midy,midz,meshType,nodesPerElement);
		
		cellIDx = floor( (midx - minx) / cellLength[0] ); 
		cellIDy = floor( (midy - miny) / cellLength[1] ); 
		cellIDz = 0; //for tri elements// floor( (midz - minz) / cellLength[2] ); 
		//cout<< " element " << element<< "cell idx " <<  cellIDx << "cell idy" << cellIDy <<endl;
		if (cellIDx < cellsNumber[0] && cellIDx>=0  && cellIDy < cellsNumber[1] && cellIDy>=0 && cellIDz < cellsNumber[2]&& cellIDz>=0)
		{
			//cout<< " element " << element<< "cell idx " <<  cellIDx << "cell idy" << cellIDy <<endl;
			// lets just use a 2d grid for the triangles
			//CellID = cellIDx * cellsNumber[1] * cellsNumber[2] + cellIDy * cellsNumber[2] + cellIDz;
			CellID = cellIDx * cellsNumber[1]  + cellIDy ;
			if (CellID>=0 && CellID < numCells)
				{
					
					//cout<< " midx " <<  midx << "midy" << midy <<endl;

					//cout<<"CellID "<<CellID<<endl;
					
					CellsVector[CellID].push_back(element);
					
					elementsFound++;
				}
		}
		//Cells.push_back(tempCell);

	}

	//cout<< " Total elements in mesh " <<numElements << "  elements in the grid " << elementsFound<<endl ;

	// go through each hinge find in which cell it is and test the hinge with all elements that are in that cell and its neighbouring cells.

	for (int fiberID= 0; fiberID < fibers.size(); fiberID++)
	{
		for (int hingeID = 0; hingeID < fibers[fiberID].numberOfHinges; hingeID++)
		{
			//cout<<"fiber " << fiberID << "  hinge " << hingeID<<endl;
			double hingeX = fibers[fiberID].hinges[hingeID].position[0];
			double hingeY = fibers[fiberID].hinges[hingeID].position[1];
			double hingeZ = fibers[fiberID].hinges[hingeID].position[2];
			//cout<< " hinge coord X " << hingeX << " hinge coord y " << hingeY <<endl;  
			// get the index of the cell where the hinge is

				cellIDx = floor( (hingeX- minx) / cellLength[0] ); 
				cellIDy = floor( (hingeY - miny) / cellLength[1] ); 
				cellIDz = floor( (hingeZ - minz) / cellLength[2] ); 
				// lets just use a 2d grid for the triangles
			//CellID = cellIDx * cellsNumber[1] * cellsNumber[2] + cellIDy * cellsNumber[2] + cellIDz;
			CellID = cellIDx * cellsNumber[1]  + cellIDy ;
				
				int neighbouringCell;
				//test with all elements in neighbouring elements

				for (int indexX  = cellIDx-1; indexX < cellIDx+2; indexX++)
				{
					
					//cout<< " index X " << indexX <<   endl; 
					if (indexX >=0 && indexX< cellsNumber[0])
					{

					
					for (int indexY  = cellIDy-1; indexY < cellIDy+2; indexY++)
						{
							//cout<< " index Y " << indexY <<"  max index Y" << cellsNumber[1] <<   endl; 
							if (indexY >=0 && indexY< cellsNumber[1])
							{
							for (int indexZ  = cellIDz-1; indexZ < cellIDz+2; indexZ++)
								{
									//cout<< " index Z " << indexZ <<"  max index Z " << cellsNumber[2] <<   endl; 
									if (indexZ >=0 && indexZ< cellsNumber[2])
										{
									//neighbouringCell = indexX * cellsNumber[1] * cellsNumber[2] + indexY * cellsNumber[2] + indexZ;
											neighbouringCell = indexX * cellsNumber[1] + indexY ; //2D
									//cout<<" test elements form cell ID " << neighbouringCell<< " hinge Cell is  "<< CellID<<endl;
									// test point with elements in this cell

									/*for (int element = 0; element < CellsVector[neighbouringCell].size(); Cells[neighbouringCell].elementID.size(); element++)
									{
										int elementID = Cells[neighbouringCell].elementID[element ];
										*/
											//cout<<"hinge Indices " << cellIDx << " "<<cellIDy<<endl;
											//cout<< " Cell Indices "<< indexX<<" " <<indexY<<endl;
										for (int element = 0; element < CellsVector[neighbouringCell].size(); element++)
									{
										int elementID = CellsVector[neighbouringCell][element ];
										//cout<< " "<< elementID<<" ";

										if(meshType =TYPE_TRI){

											bool printdebug =false;

											/*if (elementID==522)
											{
												cout<<endl<< " hinge coord X " << hingeX <<endl;  
					                          cout<< " hinge coord y " << hingeY ;  
											 cout<<" x1 " <<  velMeshPoints[  velMeshConnectivity[elementID*6+0] + numnodes * 0     ] 
											<< " x2 " <<  velMeshPoints[  velMeshConnectivity[elementID*6+1] + numnodes * 0     ] 
											<< " x3 "<<velMeshPoints[  velMeshConnectivity[elementID*6+2] + numnodes * 0     ]
											<< " y1 " << velMeshPoints[  velMeshConnectivity[elementID*6+0] + numnodes * 1     ] 
											<< " y2 " << velMeshPoints[  velMeshConnectivity[elementID*6+1] + numnodes * 1     ] 
											<< " y3 "<<velMeshPoints[  velMeshConnectivity[elementID*6+2] + numnodes * 1     ] 
											<<endl;
											cout<< " element nodes "<<  velMeshConnectivity[elementID*6+0] << " " << velMeshConnectivity[elementID*6+1]<< " "<< velMeshConnectivity[elementID*6+2]<<endl;
											printdebug =true;
											}
											*/
											/*cout<< " hinge coord X " << hingeX <<endl;  
					                          cout<< " hinge coord y " << hingeY ;  
											 cout<<" x1 " <<  velMeshPoints[  velMeshConnectivity[elementID*6+0] + numnodes * 0     ] 
											<< " x2 " <<  velMeshPoints[  velMeshConnectivity[elementID*6+1] + numnodes * 0     ] 
											<< " x3 "<<velMeshPoints[  velMeshConnectivity[elementID*6+2] + numnodes * 0     ]
											<< " y1 " << velMeshPoints[  velMeshConnectivity[elementID*6+0] + numnodes * 1     ] 
											<< " y2 " << velMeshPoints[  velMeshConnectivity[elementID*6+1] + numnodes * 1     ] 
											<< " y3 "<<velMeshPoints[  velMeshConnectivity[elementID*6+2] + numnodes * 1     ] 
											<<endl;*/

											//cout<< " element nodes "<<  velMeshConnectivity[elementID*6+0] << " " << velMeshConnectivity[elementID*6+1]<< " "<< velMeshConnectivity[elementID*6+2]<<endl;
											if(testPointandTri(hingeX,hingeY,0,
												velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+0] + numnodes * 0     ], //first node(0) position x(0) 
												velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+0] + numnodes * 1     ], //first node(0) position y(1) 
												0, //first node(0) position z(2) 
												velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+1] + numnodes * 0     ], //second node(1) position y(0) 
												velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+1] + numnodes * 1     ], //second node(1) position y(1) 
												0, //second node(1) position z(2) 
												velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+2] + numnodes * 0     ], //third node(2) position y(0) 
												velMeshPoints[  velMeshConnectivity[elementID*nodesPerElement+2] + numnodes * 1     ], //third node(2) position y(1) 
												0 //third node(21) position z(2) 
												,printdebug))
												{  fibers[fiberID].hinges[hingeID].elementID = elementID;  
												//cout<<"fiber " << fiberID << " hinge " << hingeID  << " is in element " << fibers[fiberID].hinges[hingeID].elementID <<endl  ;
												
												}else{//cout<<"hinge is out of element " <<elementID<<endl;
											}
 										}// if

									}// for element

										//cout<<endl;

									}else{//cout<<"hinge is out of z domain"<<endl;
									}// if indexZ
								}// for indexZ
							}else{//cout<<"hinge is out of y domain"<<endl;
							}//if indexY
						}//for indexY
					}else{//cout<<"hinge is out of x domain"<<endl;
					}//ifindexX
				}//for indexX



		}


	}

	




}

void findBoundaryEdges(int numElements, int numFacesPerElement, vector<int>& elementNeighbors, vector<int > &boundaryFaces){

	for (int elementID = 0; elementID < numElements; elementID++)
	{
		boundaryFaces[elementID]=-1;
		for (int faceID = 0; faceID < numFacesPerElement; faceID++)
		{
			if ( elementNeighbors[ elementID*numFacesPerElement + faceID] == -1 )
			{
				boundaryFaces[elementID] = faceID;
			}


		}



	}



}
