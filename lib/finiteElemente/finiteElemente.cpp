#include "finiteElemente.h"


FiniteElemente::FiniteElemente(double len_, double width_, int maxNumberOfElments){
    len=len_;
    width=width_;

    initMesh(maxNumberOfElments);

    // 5. --> see see mfem ex1.cpp in mfem lib
    fec = new H1_FECollection(order, dim);
    fespace = new FiniteElementSpace(mesh, fec);

    // 8.  --> see see mfem ex1.cpp in mfem lib
    solutionVector= new GridFunction(fespace); // changed to pointer
    *solutionVector = 0;
}


void FiniteElemente::initMesh(int maxNumberOfElments){
    double area=len*width;
    numberVerticesX=int(len/std::sqrt(area/maxNumberOfElments))+1;
    numberVerticesY=int(width/std::sqrt(area/maxNumberOfElments))+1;

    // std::cout<<"numberVerticesX "<<numberVerticesX<<" numberVerticesY "<<numberVerticesY<<" pro "<<numberVerticesX*numberVerticesY<<std::endl;


    int nv=numberVerticesX*numberVerticesY;
    int ne=(numberVerticesX-1)*(numberVerticesY-1);
    int nb=2*(numberVerticesX-1)+2*(numberVerticesY-1);

    mesh=new Mesh(dim, nv, ne, nb, sdim);


    // add elements
    vertexIndexMap = new int*[numberVerticesX];
    for(int i=0;i<numberVerticesX;i++){
        vertexIndexMap[i]= new int[numberVerticesY];
        for(int j=0;j<numberVerticesY;j++){
            mesh->AddVertex(Vertex(len*i/(numberVerticesX-1),width*j/(numberVerticesY-1))()); 
            vertexIndexMap[i][j]=numberVerticesY*i+j;
        }
    }
    
    // add elements
    int quadIndx[4];
    for(int i=0;i<numberVerticesX-1;i++){
        for(int j=0;j<numberVerticesY-1;j++){
            quadIndx[0]=vertexIndexMap[i  ][j  ];
            quadIndx[1]=vertexIndexMap[i+1][j  ];
            quadIndx[2]=vertexIndexMap[i+1][j+1];
            quadIndx[3]=vertexIndexMap[i  ][j+1];

            // std::cout<<vertexIndexMap[i][j]<<vertexIndexMap[i+1][j]<<vertexIndexMap[i+1][j+1]<<vertexIndexMap[i][j+1]<<std::endl;
            mesh->AddQuad(quadIndx);
        }
    }

    // add boundaries
    int bdrIndx[2];
    for(int i=0;i<numberVerticesX-1;i++){
        bdrIndx[0]=vertexIndexMap[i+1][0];
        bdrIndx[1]=vertexIndexMap[i  ][0];
        mesh->AddBdrSegment(bdrIndx);
        bdrIndx[0]=vertexIndexMap[i  ][numberVerticesY-1];
        bdrIndx[1]=vertexIndexMap[i+1][numberVerticesY-1];
        mesh->AddBdrSegment(bdrIndx);
    }
    for(int j=0;j<numberVerticesY-1;j++){
        bdrIndx[0]=vertexIndexMap[0][j  ];
        bdrIndx[1]=vertexIndexMap[0][j+1];
        mesh->AddBdrSegment(bdrIndx);
        bdrIndx[0]=vertexIndexMap[numberVerticesX-1][j+1];
        bdrIndx[1]=vertexIndexMap[numberVerticesX-1][j  ];
        mesh->AddBdrSegment(bdrIndx);
    }
    mesh->FinalizeQuadMesh(0, 0, true);

}

void FiniteElemente::setElectrode(double begin, double end,int edge,double voltage){
    // set electrode in given range on edge given by:
    // 0 --> x=0
    // 1 --> x=len
    // 2 --> y=0
    // 3 --> y=width


    //validity checks
    if(edge > 3 || edge < 0){
        throw std::invalid_argument("invalid edge given");
    }
    if((edge == 0 || edge == 1) && (end > width ||  begin < 0)){
        throw std::length_error("electrode out of range");
    }
    else if((edge == 2 || edge == 3) && (end > len ||  begin < 0)){
        throw std::length_error("electrode out of range");
    }


    //get vertex indices inside range
    std::vector<int> vertexIndices;
    switch (edge){
        case 0:
            for(int j=0;j<numberVerticesY;j++){
                if(mesh->GetVertex(vertexIndexMap[0][j])[1]>=begin and mesh->GetVertex(vertexIndexMap[0][j])[1]<=end){
                    vertexIndices.push_back(j);
                }
            }
            break;
        case 1:
            for(int j=0;j<numberVerticesY;j++){
                if(mesh->GetVertex(vertexIndexMap[numberVerticesX-1][j])[1]>=begin and mesh->GetVertex(vertexIndexMap[numberVerticesX-1][j])[1]<=end){
                    vertexIndices.push_back(vertexIndexMap[numberVerticesX-1][j]);
                }
            }
            break;
        case 2:
            for(int i=0;i<numberVerticesX;i++){
                if(mesh->GetVertex(vertexIndexMap[i][0])[0]>=begin and mesh->GetVertex(vertexIndexMap[i][0])[0]<=end){
                    vertexIndices.push_back(vertexIndexMap[i][0]);
                }
            }
            break;
        case 3:
            for(int i=0;i<numberVerticesX;i++){
                if(mesh->GetVertex(vertexIndexMap[i][numberVerticesY-1])[0]>=begin and mesh->GetVertex(vertexIndexMap[i][numberVerticesY-1])[0]<=end){
                    vertexIndices.push_back(vertexIndexMap[i][numberVerticesY-1]);
                }
            }
            break;
    }



    // set initial condition in solution vector
    for(auto index:vertexIndices){
        (*solutionVector)[index]=voltage;
    }



    // set boundary attribute to make (only) it mandatory
    Array<int> bdrVertices = {0,0};
    for(int i=0;i<2*(numberVerticesX-1)+2*(numberVerticesY-1);i++){ //loop over all boundary elements
        mesh->GetBdrElementVertices(i,bdrVertices);
        for(auto index1:vertexIndices){
            if (bdrVertices[0] == index1){
                for(auto index2:vertexIndices){
                    if (bdrVertices[1] == index2){
                        mesh->GetBdrElement(i)->SetAttribute(2);
                        break;
                    }
                }
            }
        }    
    }


    
    // move first and last vertex to exactly fit electrode boundaries
    if (vertexIndices.size()<2){
        throw std::length_error("could not find 2 vertices inside electrode range, refine grid or enlarge electrode");
        }
    else{
        switch (edge){
            case 0:
            case 1:
                mesh->GetVertex(vertexIndices.front())[1]=begin;
                mesh->GetVertex(vertexIndices.back())[1]=end;
                break;
            case 2:
            case 3:
                mesh->GetVertex(vertexIndices.front())[0]=begin;
                mesh->GetVertex(vertexIndices.back())[0]=end;
                break;
                
            }
        }
}

void FiniteElemente::run(){
    // see mfem ex1.cpp in mfem lib


    // 2.--> see see mfem ex1.cpp in mfem lib
   Device device("cpu");
//    device.Print();



   //6.--> see see mfem ex1.cpp in mfem lib
   Array<int> ess_tdof_list;
   if (mesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh->bdr_attributes.Max());
    //   ess_bdr = 0;
      ess_bdr[1] = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
   //    the basis functions in the finite element fespace.
   LinearForm *b = new LinearForm(fespace);
   ConstantCoefficient one(1.0);
   ConstantCoefficient zero(0.0);
   b->AddDomainIntegrator(new DomainLFIntegrator(zero));
   b->Assemble();

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm *a = new BilinearForm(fespace);
   a->AddDomainIntegrator(new DiffusionIntegrator(one));

   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.
   a->Assemble();

   OperatorPtr A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, *solutionVector, *b, A, X, B);


   // 11. Solve the linear system A X = B.
    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, X, 0, 1000, 1e-12, 0.0);

   // 12. Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(X, *b, *solutionVector);


   // -------- tetsing only -------
   // 13. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
   ofstream mesh_ofs("refined.mesh");
   mesh_ofs.precision(8);
   mesh->Print(mesh_ofs);
   ofstream sol_ofs("sol.gf");
   sol_ofs.precision(8);
   solutionVector->Save(sol_ofs);

}

double FiniteElemente::getPotential(double x, double y){
    if(x>len) throw std::invalid_argument("x pos of getPotential out of range");
    if(y>len) throw std::invalid_argument("y pos of getPotential out of range");
    return (*solutionVector)[vertexIndexMap[int(x/len*(numberVerticesX-1)+0.5)][int(y/width*(numberVerticesY-1)+0.5)]];
}


FiniteElemente::~FiniteElemente(){
    delete fespace;
    delete fec;
    delete mesh;
    delete vertexIndexMap;
}
