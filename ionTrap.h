#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "vector.h"

#include "TFile.h"
#include "TTree.h"


class ionTrap {
    
public:
    double Vac_min = 0;
    double Vac_max = 500;
    double Vac_step = 5;
    double Vac;
    
    double Vdc_min = -150;
    double Vdc_max = 150;
    double Vdc_step = 3;
    double Vdc;
    
    double w_min = 1e5;
    double w_max = 9e5;
    double w_step =1e5;
    double w;
    
    vector Fg;
    
    // For the Drag force (in kg/s)
    double  drag = 0;
    //    double drag = 3.47e-16;
    
    // For the Brownian kick (in J)
    double kT = 0;
    //double kT = 1.2e-21;
    
    // For the eletric force
    double  q = -1.6e-19;
    double  z0 = 3e-3;
    
    // Particle mass (in kg)
    double  m = 2.28e-25;
    
    // Time step (in s)
    double  dt = 1e-8;
    double  t = 0;
    // Brownian kick
    double  kickvar = 2*kT*drag/dt;
    
    double  a_z = -16*Vdc*q/(m*z0*z0*w*w);
    double  q_z = 8*Vac*q/(m*z0*z0*w*w);
    
    // Position and velocity (in m and m/s)
    vector r;
    vector v;
    double v_terminal = mag(Fg)/drag;
    double sumVelocity;
    double sumForce;
    
    double timeCounter = 0;
    double position;
    double positionZ;
    int counter = 0;
    
    std::string rootOutputFile = "trapOutput_vacuum_e7.root";
    
    //------------------------------------------//
    // Class methods
    
    ionTrap() {
        r.reassign(0.001, 0.0005, 0.001);
        v.reassign(-0.5, 0.5, -0.1);
    }
    
    ~ionTrap() {}
    
    // Randomly generating the sign of the Brownian force term
    double  kick_sign() {
        
        int random = rand() % 100 + 1;
        
        if (random % 2 == 0) {
            return sqrt(kickvar);
        }
        else {
            return -sqrt(kickvar);
        }
    }
    
    // Calculating the force vector the particle is under
    vector F(vector r, vector v, double  t){
        
        //The eletric term
        static vector Fe = q * vector((Vdc - Vac * cos (w * t)) / (z0*z0) * ( 2 * r.x),
                                      (Vdc - Vac * cos (w * t)) / (z0*z0) * ( 2 * r.y),
                                      (Vdc - Vac * cos (w * t)) / (z0*z0) * (-1 * r.z));
        //The drag term
        static vector Fd = -1 * v * drag;
        
        //The Brownian term
        static vector Fb = vector((kick_sign() * kickvar), (kick_sign() * kickvar), (kick_sign() * kickvar));
        
        static vector Fg = vector(0,-9.8*m, 0);
        
        return Fe + Fd + Fb + Fg;
    }
    
    
    void ExecuteLeapFrog(){
        
         r += v * dt/2;
         v += F(r,v,t)/m * dt/2;
         v += F(r,v,t + dt/2)/m * dt/2;
         r += v * dt/2;
         t += dt;
    }
    
    
    void BuildRootFile() {
        // Defining vectors to be used in the root file
        std::vector<double> vector_X;
        std::vector<double> vector_Y;
        std::vector<double> vector_Z;
        std::vector<double> vector_t;
        std::vector<double> vector_vX;
        std::vector<double> vector_vY;
        std::vector<double> vector_vZ;
        std::vector<double> vector_Vac;
        std::vector<double> vector_Vdc;
        std::vector<double> vector_frequency;
        std::vector<double> vector_normR;
        std::vector<double> vector_normZ;
        std::vector<double> vector_avgV;
        std::vector<int> vector_counter;
        std::vector<double> vector_avgF;
        std::vector<double> vector_finalF;
        
        // Creating the root file
        TFile rootOutput((const char*) rootOutputFile.c_str(), "RECREATE");
        
        TTree inputTree("inputs", "inputs");
        
        inputTree.Branch("Vac_min", &Vac_min, "Vac_min/D");
        inputTree.Branch("Vac_max", &Vac_max, "Vac_max/D");
        inputTree.Branch("Vac_step", &Vac_step, "Vac_step/D");
        inputTree.Branch("Vdc_min", &Vdc_min, "Vdc_min/D");
        inputTree.Branch("Vdc_max", &Vdc_max, "Vdc_max/D");
        inputTree.Branch("Vdc_step", &Vdc_step, "Vdc_step/D");
        inputTree.Branch("w_min", &w_min, "w_min/D");
        inputTree.Branch("w_max", &w_max, "w_max/D");
        inputTree.Branch("w_step", &w_step, "w_step/D");
        inputTree.Branch("brownianKick", &kT, "brownianKick/D");
        inputTree.Branch("drag", &drag, "drag/D");
        inputTree.Branch("q", &q, "q/D");
        inputTree.Branch("z0", &z0, "z0/D");
        inputTree.Branch("v_terminal", &v_terminal, "v_terminal/D");
        
        TTree trapTree("trap", "trap");
        
        trapTree.Branch("X", &vector_X);
        trapTree.Branch("Y", &vector_Y);
        trapTree.Branch("Z", &vector_Z);
        trapTree.Branch("Tempo", &vector_t);
        trapTree.Branch("norm_R", &vector_normR);
        trapTree.Branch("norm_Z", &vector_normZ);
        trapTree.Branch("vX", &vector_vX);
        trapTree.Branch("vY", &vector_vY);
        trapTree.Branch("vZ", &vector_vZ);
        trapTree.Branch("Vac", &vector_Vac);
        trapTree.Branch("Vdc", &vector_Vdc);
        trapTree.Branch("w", &vector_frequency);
        
        trapTree.Branch("avgV", &vector_avgV);
        trapTree.Branch("avgForce", &vector_avgF);
        trapTree.Branch("finalForce", &vector_finalF);
        trapTree.Branch("#_cycles_before_ending", &vector_counter);
        
        
        for (Vac = Vac_min; Vac < Vac_max; Vac += Vac_step) {
            
            for (Vdc = Vdc_min; Vdc < Vdc_max; Vdc += Vdc_step){
                
                for (w = w_min; w < w_max; w += w_step) {
                    
                    // Calculating kickvar
                    kickvar = 2 * kT * drag/dt;
                    
                    while (t < (10000/w)) {
                        
                        ExecuteLeapFrog();
                        sumVelocity += mag(v);
                        sumForce += mag(F(r,v,t));
                        
                        counter++;
                        
                        if(mag(r)> 1.42*z0){
                            continue;
                        }
                    }
                    
                    double normZ = sqrt(r.z*r.z);
                    double avgV = sumVelocity/counter;
                    double avgF = sumForce/counter;
                    
                    
                    vector_X.push_back(r.x);
                    vector_Y.push_back(r.y);
                    vector_Z.push_back(r.z);
                    vector_t.push_back(t);
                    
                    vector_vX.push_back(v.x);
                    vector_vY.push_back(v.y);
                    vector_vZ.push_back(v.z);
                    
                    vector_Vac.push_back(Vac);
                    vector_Vdc.push_back(Vdc);
                    vector_frequency.push_back(w);
                    
                    vector_normR.push_back(norm(r));
                    vector_normZ.push_back(normZ);
                    
                    vector_avgV.push_back(avgV);
                    vector_avgF.push_back(avgF);
                    vector_finalF.push_back(mag(F(r,v,t)));
                    vector_counter.push_back(counter);
                    
                    r.reassign(0.001, 0.0005, 0.001);
                    v.reassign(-0.5, 0.5, -0.1);
                    
                    t = 0;
                    counter = 0;
                    
                }
            }
        }
        // For the root file
        inputTree.Fill();
        inputTree.Write();
        trapTree.Fill();
        trapTree.Write();
        rootOutput.Close();
    }
    
    void PrintForAnime() {
        double time_factor = 1e-5; // fraction of realtime to run at
        int frameskip = 1./(60 * dt) * time_factor;
        printf("!Beginning simulation. Frameskip = %d\n",frameskip);
        printf("!Parameters:\n");
        printf("!  Vac = %e \t Vdc = %e \t w = %e \t dt = %e \t drag = %e \t kT = %e\n",Vac,Vdc,w,dt,drag,kT);
        printf("!  a_z = %e \t q_z = %e\n",a_z,q_z);
        printf("!Sanity checks: dimensionless period %e | drag decay constant %e\n",dt*w,drag/m * dt);
        int steps;
        while (true)  {
            ExecuteLeapFrog();
            steps++;
            if (steps % frameskip == 0)
            {
                printf("ct3 0 %e %e %e %e\n", r.x, r.y, r.z, 1e-4);
                printf("F\n");
            }
        }
    }
};




