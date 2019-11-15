#include "finiteElemente.h"


FiniteElemente::FiniteElemente(double len_, double width_, int maxNumberOfElments){
    len=len_;
    width=width_;


    double area=len*width;
    numberVerticesX=int(len/std::sqrt(area/maxNumberOfElments))+1;
    numberVerticesY=int(width/std::sqrt(area/maxNumberOfElments))+1;

    std::cout<<"numberVerticesX "<<numberVerticesX<<" numberVerticesY "<<numberVerticesY<<" pro "<<numberVerticesX*numberVerticesY<<std::endl;


    int nv=numberVerticesX*numberVerticesY;
    int ne=(numberVerticesX-1)*(numberVerticesY-1);
    int nb=2*(numberVerticesX-1)+2*(numberVerticesY-1);

    mesh=new Mesh(dim, nv, ne, nb, sdim);


    for(int i=0;i<numberVerticesX;i++){
        for(int j=0;j<numberVerticesY;j++){
            mesh->AddVertex(Vertex(len*i/(numberVerticesX-1),width*j/(numberVerticesY-1))()); 
            // std::cout<<"v_x: "<<len*i/(numberVerticesX-1)<<" v_y: "<<width*j/(numberVerticesY-1)<<std::endl;
        }
    }
    
    int quadIndx[4];
    for(int i=0;i<numberVerticesX-1;i++){
        for(int j=0;j<numberVerticesY-1;j++){
            quadIndx[0]=i*    numberVerticesY+j;
            quadIndx[1]=(i+1)*numberVerticesY+j;
            quadIndx[2]=(i+1)*numberVerticesY+(j+1);
            quadIndx[3]=i*    numberVerticesY+(j+1);

            // std::cout<<quadIndx[0]<<quadIndx[1]<<quadIndx[2]<<quadIndx[3]<<std::endl;
            mesh->AddQuad(quadIndx);
        }
    }

    int bdrIndx[2];
    for(int i=0;i<numberVerticesX-1;i++){
        bdrIndx[0]=(i+1)*numberVerticesY;
        bdrIndx[1]=i*numberVerticesY;
        mesh->AddBdrSegment(bdrIndx);
        bdrIndx[0]=i    *numberVerticesY+(numberVerticesY-1);
        bdrIndx[1]=(i+1)*numberVerticesY+(numberVerticesY-1);
        mesh->AddBdrSegment(bdrIndx,2);
    }
    for(int j=0;j<numberVerticesY-1;j++){
        bdrIndx[0]=j;
        bdrIndx[1]=j+1;
        mesh->AddBdrSegment(bdrIndx,3);
        bdrIndx[0]=j+1 +numberVerticesY*(numberVerticesX-1);
        bdrIndx[1]=j   +numberVerticesY*(numberVerticesX-1);
        mesh->AddBdrSegment(bdrIndx,4);
    }
    mesh->FinalizeQuadMesh(0, 0, true);

}

void FiniteElemente::setElectrode(double begin, double end,int edge,double voltage){
    // set electrode in given range on edge given by:
    // 0 --> x=0
    // 1 --> x=len
    // 2 --> y=0
    // 3 --> y=width

    //get vertex indices inside range
    std::vector<int> vertexIndices;
    switch (edge){
        case 0:
            for(int j=0;j<numberVerticesY;j++){
                if(mesh->GetVertex(j)[1]>=begin and mesh->GetVertex(j)[1]<=end){
                    vertexIndices.push_back(j);
                }
            }
            break;
        case 1:
            for(int j=0;j<numberVerticesY;j++){
                if(mesh->GetVertex(j)[1]>=begin and mesh->GetVertex(j)[1]<=end){
                    vertexIndices.push_back(j);
                }
            }
            break;
    }


    std::cout<<"found: ";
    for(auto a:vertexIndices){
        std::cout<<a<<" ";
    }
    std::cout<<std::endl;


    
    // move first and last vertex
    if (vertexIndices.size()<2){
        throw std::invalid_argument("could not find 2 vertices inside electrode range, refine grid or enlarge electrode");
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
    // see ex1.cpp in mfem lib


    // 2.
   Device device("cpu");
   device.Print();


    // 5.
   FiniteElementCollection *fec;
   int order =1;
   fec = new H1_FECollection(order, dim);

   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);

   //6.
   Array<int> ess_tdof_list;
   if (mesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }



   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
   //    the basis functions in the finite element fespace.
   LinearForm *b = new LinearForm(fespace);
   ConstantCoefficient one(1.0);
   ConstantCoefficient zero(0.0);
   b->AddDomainIntegrator(new DomainLFIntegrator(one));
   b->Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = 0;
   for(int i=0;i<20;i++){
       x[i+10]=10;
   }
   
   std::cout<<"x size"<<x.Size()<<std::endl;

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
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);


   // 11. Solve the linear system A X = B.
    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, X, 1, 200, 1e-12, 0.0);

   // 12. Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(X, *b, x);


   // -------- tetsing only -------
   // 13. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
   ofstream mesh_ofs("refined.mesh");
   mesh_ofs.precision(8);
   mesh->Print(mesh_ofs);
   ofstream sol_ofs("sol.gf");
   sol_ofs.precision(8);
   x.Save(sol_ofs);


    // -------- tetsing only -------
    // 14. Send the solution by socket to a GLVis server.
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sol_sock(vishost, visport);
    sol_sock.precision(8);
    sol_sock << "solution\n" << *mesh << x << flush;


    delete fespace;
    delete fec;
}



void FiniteElemente::printMesh(std::string fileName){

    ofstream meshfile;

    meshfile.open(fileName,ios::trunc);
    mesh->Print(meshfile);
    meshfile.close();

}

FiniteElemente::~FiniteElemente(){
   delete mesh;

}
