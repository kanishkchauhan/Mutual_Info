/**
 *        File: firingVScoupling.cc
 *      Author: Kanishk Chauhan
 *        Date: June 8, 2019
 * Description: to determine the dependence of firing rate and the coefficient of variance 
 *              on the coupling strength and the number of nodes in a Hodgkin-
 *              huxley model.
 *               
 */


# include <iostream>
# include <cmath>
# include <vector>
# include <cstdlib> //needed for random numbers
#include <fstream> // needed to write the threshold value data to a text file

using namespace std;

// the following funtions will calculate changes in V,m,a nd h for each iteration
double hhnode_fV(double v, double v_central, double m, double h, double Iext, double kappa);
double hhnode_fV_central(double v, double v_peri_nodes[], double m, double h, double kappa, int N);
double hhnode_fm(double v, double m);
double hhnode_fh(double v, double h);

const double tmax = 200; //max time in ms (equal to 10 s)
const double dt = 0.001; //time step
const double tu = 100; //relaxation time

int main ()
{
    // N is the number of nodes including the central one, the Nth node being the central one
    // kappa is the coupling constant

    //following chunk of codes is needed to increment the value of kappa logarithmically 
    double kappa_min = 0.4;
    double kappa_max = 1000;
    double number_of_points = 30;
    double alpha = 1/(number_of_points - 1) * log10(kappa_max/kappa_min);
    // this 'alpha' will be used to obtain the next value of kappa

    ofstream file; //creating the file to store data
    file.open("threshold_VS_coupling_data.csv");
    file << "Coupling,F-R central,CV central,F-R peripheral,CV peripheral" << endl;

    for (int N = 2; N <= 10; N++) 
    {
        file << endl;
        file << "For N = " << N << endl;
        double v0[N], m0[N], h0[N]; // arrays for initial conditions
        double variation_amp_v0[N], variation_amp_m0[N], variation_amp_h0[N];
        double Iext = 20.0; // starting value of external current- 
                                // this will be increased in small steps in the following while loop 
        
        double kappa = kappa_min; // initialization for while loop 
        int counter = 1; //to obtain next value of kappa (logarithmically) this is used as 'i'

        while (kappa < kappa_max)
        {
            
            double spike_time;
            // the following vector v is to initiate the 2D vector spike_mat to store the spike times of all N nodes
            vector<double> v;
            vector<vector<double> > spike_mat(N); // creates an N by 1 2D matrix/vector - we can add elements and values later

            counter ++; // this is for incrementing the value of kappa
            
            for (int i =0; i <= N-1; i++)
            {
                // we begin with generating different initial values of V,m,and h for the nodes
                v0[i] = -80; m0[i] = 0.2; h0[i] = 0.5;
                variation_amp_v0[i] = rand()%20-10 ; // small variations in values of v for nodes 
                variation_amp_m0[i] = (rand()%10-5) * 0.01;//small variations in values of m for nodes
                variation_amp_h0[i] = (rand()%10-5) * 0.01;//small variations in values of h for nodes
                v0[i] = v0[i] + variation_amp_v0[i];
                m0[i] = m0[i] + variation_amp_m0[i];
                h0[i] = h0[i] + variation_amp_h0[i];
            } // creates different initial conditions for all nodes by adding small variations
              // for each value of N, we will examine all possible values of kappa
            
            // we start with transition to steady state
            double nu = tu/dt + 1; 
            for (double j = 1; j <= nu; j += 1) // this updates the values of V,m,and h
            {   
                double v_central = v0[N-1];
                for (int k = 0; k <= N-1; k++)
                {
                    double v = v0[k], m = m0[k], h = h0[k]; // assigns initial values of V,m, and h corresponding to the kth node to new variables to pass through funtions                        double v_central = v0[N-1];
                    double fV; 
                    if (k == N-1)
                    { 
                        double v_peri_nodes[N-1];
                        for (int l = 0; l < N-1; l++)
                        {
                            v_peri_nodes[l] = v0[l];
                        }
                        fV = hhnode_fV_central(v,v_peri_nodes,m,h,kappa,N); 
                    } // calculates the change in value of V for central node
                    else { fV = hhnode_fV(v,v_central,m,h,Iext,kappa); } // calculates the change in value of V
                    double fm = hhnode_fm(v,m); // calculates the change in value of m
                    double fh = hhnode_fh(v,h); // calculates the change in value of h

                    v0[k] = v0[k] + (dt * fV); // Euler steps
                    m0[k] = m0[k] + (dt * fm);
                    h0[k] = h0[k] + (dt * fh);
                } // transition to steady ends here for kth node
            } // transition to steady ends here for all N nodes

            double nt = tmax/dt; 
            double v_new[N]; double m_new[N]; double h_new[N]; //needed for updating the values 
                                                                //of V,m,and h                                                                
            int iflag = 1;
            for (double jj = 1; jj <= nt; jj++) // this updates the values of V,m,and h for each node one by one
            {
                double v_central = v0[N-1];
               
                for (int kk = 0; kk <= N-1; kk++)
                {
                    double v = v0[kk], m = m0[kk], h = h0[kk]; // assigns initial values for kth node to new variables to pass through funtions
                    double fV;
                
                    if (kk == N-1) 
                    {
                        double v_peri_nodes[N-1];
                        for (int l = 0; l < N-1; l++)
                        {
                            v_peri_nodes[l] = v0[l];
                        }
                        fV = hhnode_fV_central(v,v_peri_nodes,m,h,kappa,N);
                    } // calculates the change in value of V for central node
                    else { fV = hhnode_fV(v,v_central,m,h,Iext,kappa); } // calculates the change in value of V
                    double fm = hhnode_fm(v,m); // calculates the change in value of m
                    double fh = hhnode_fh(v,h); // calculates the change in value of h

                    v_new[kk] = v0[kk] + (dt * fV); // Euler steps
                    m_new[kk] = m0[kk] + (dt * fm);
                    h_new[kk] = h0[kk] + (dt * fh);
                   
                    // discriminating spikes for all node
                    if (v0[kk] > -20 && v_new[kk] < -20)
                    {
                        iflag = 1;
                    }
                    if (v0[kk] < 20 && v_new[kk] > 20 && iflag == 1)
                    {
                        spike_time = (jj - 1)*dt; // this stores spike time for kkth node
                        spike_mat[kk].push_back(spike_time);
                        iflag = 0;
                    }
                    //updating the values of V,m, and h
                    v0[kk] = v_new[kk];
                    m0[kk] = m_new[kk];
                    h0[kk] = h_new[kk];                              }                     
                }    
            }// spike_mat now has all spike time values for all nodes 

            // calculations for inter-spike interval and firing rate
            // the following vectors will store the isi for one of the peripheral nodes and the central node
            vector<double> isi_peri(1),isi_central(1);
            double sum1 = 0, sum2 = 0;
            double meanISI_peri, meanISI_central, firing_rate_peri, firing_rate_central, CV_peri, CV_central;
           // int number_spikes_peri = spike_mat[0].size();
           //int number_spikes_central = spike_mat[N-1].size();
            for (size_t val = 0; val < spike_mat[1].size()-1; val++)
            {
                isi_peri[val] = spike_mat[1][val+1] - spike_mat[1][val];
                sum1 += isi_peri[val];
            }
            meanISI_peri = sum1 / (isi_peri.size());
            firing_rate_peri = 1000/meanISI_peri;
            CV_peri = 0;

            for (size_t vall = 0; vall < spike_mat[N-1].size()-1; vall++)
            {
                isi_central[vall] = spike_mat[N-1][vall+1] - spike_mat[N-1][vall];
                sum2 += isi_central[vall];
            }
            meanISI_central = sum2 / (isi_central.size());
            firing_rate_central = 1000/meanISI_central;
            CV_central = 0;
            //cout << "The required threshold current is: " << Iext-0.1 << " for kappa "<< kappa << endl; // because last iteration of the previous loop adds an extra 0.1 to Iext
            file << kappa << "," << firing_rate_central << "," << CV_central << "," << firing_rate_peri << "," << CV_peri << endl;
            double expo = alpha * (counter - 1);// for updating the value of kappa
            kappa = kappa_min * pow(10,expo);
        }// firing rates, etc have been calculated and written in a file            
    }
    file.close();
    return 0;
}


double hhnode_fm(double v, double m)
{
    double v1 = v + 20.4; double v2 = v + 25.7;
    double am = (1.314 * v1) / (1.0 - exp(-v1/10.3));
    double bm = -(0.0608 * v2) / (1.0 - exp(v2/9.16));
    double fm = am * (1.0 - m) - bm * m;
    return fm;
}

double hhnode_fh(double v, double h)
{
    double v3 = v + 114.0;
    double ah = -(0.068 * v3) / (1.0- exp(v3/11.0));
    double bh = 2.52 / (1.0 + exp(-(v+31.8)/13.4));
    double fh = ah * (1.0 - h) - bh * h;
    return fh;
}

double hhnode_fV(double v,double v_central, double m, double h, double Iext, double kappa)
{
    double I_ion = 20.0 * (v + 80.0) + m*m*m*h * 1100.0 * (v - 50.0);
    double fV = (-I_ion + Iext + kappa * (v_central - v))/2;
    return fV;
}

double hhnode_fV_central(double v, double v_peri_nodes[], double m, double h, double kappa, int N)
{
    double I_ion = 20.0 * (v + 80.0) + m*m*m*h * 1100.0 * (v - 50.0);
    double sum = 0;
    for (int i = 0; i < N-1; i++)
    {
        sum += v_peri_nodes[i];
    }
    double fV = (-I_ion - kappa * ((N-1)*v - sum))/2;
    return fV;
}
