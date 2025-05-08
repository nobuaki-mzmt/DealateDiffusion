// Data analysis on tired termite project

/*------------------------------------------------------------------------------
 onesim.cpp, 
 a additional code for simulations.R, 
 includes functions for calculations
 ------------------------------------------------------------------------------*/

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// Random sampling
#include <random>
#include <time.h>
#include <cmath>

// [[Rcpp::export]]
List one_simulation(double all_L,    // field size
                             double Light_L,  // light trapped area size
                             double detection,
                             int end_time,
                             int num_ind,
                             List solo_param
                             ) {
  int i_ind, i_time, j_ind;
  double dx, dy, distance;
  double u, sign_u, abs_u, laplace_value; // for laplace random number
  
  double acc_inter  = solo_param["acc_inter"];
  double acc_slope  = solo_param["acc_slope"];
  double acc_scale  = solo_param["acc_scale"];
  double turn_scale = solo_param["turn_scale"];
  
  // Initialization
  NumericMatrix population(num_ind, 6); // Columns: X, Y, Angle, Speed, Sex, Encounter
  for(i_ind = 0; i_ind < num_ind; i_ind++) {
    population(i_ind, 0) = R::runif(0, Light_L);           // X
    population(i_ind, 1) = R::runif(0, Light_L);           // Y
    population(i_ind, 2) = R::runif(0, 2 * M_PI);          // Angle
    population(i_ind, 3) = R::runif(0, 8);                 // Speed
    population(i_ind, 4) = (i_ind < num_ind / 2) ? 0 : 1;  // Sex: 0 = female, 1 = male
    population(i_ind, 5) = 0;                              // Encounter: 0 = no encounter
  }  
  
  NumericVector num_encounter(end_time);
    
  // Simulation
  for(i_time = 0; i_time < end_time; i_time++) {
    
    // Move individuals
    for(i_ind = 0; i_ind < num_ind; i_ind++) {
      if(population(i_ind, 5) == 1){ continue; }
      // Update position based on speed and angle
      population(i_ind, 0) += population(i_ind, 3) * cos(population(i_ind, 2)); // Update X
      population(i_ind, 1) += population(i_ind, 3) * sin(population(i_ind, 2)); // Update Y
      
      // Ensure individuals stay within the boundaries using wrapping
      if (population(i_ind, 0) > all_L) population(i_ind, 0) -= all_L;
      if (population(i_ind, 1) > all_L) population(i_ind, 1) -= all_L;
      if (population(i_ind, 0) < 0) population(i_ind, 0) += all_L;
      if (population(i_ind, 1) < 0) population(i_ind, 1) += all_L;
    }
    
    // Check for encounters (male and female within detection range)
    for(i_ind = 0; i_ind < num_ind / 2; i_ind++) {  // Loop over females
      if(population(i_ind, 5) == 1){ continue; }
      for (j_ind = num_ind / 2; j_ind < num_ind; j_ind++) {  // Loop over males
        if(population(j_ind, 5) == 1){ continue; }
        dx = population(i_ind, 0) - population(j_ind, 0);
        dy = population(i_ind, 1) - population(j_ind, 1);
        if (abs(dx) > all_L / 2) { dx = all_L - abs(dx); }
        if (abs(dy) > all_L / 2) { dy = all_L - abs(dy); }
        distance = sqrt(dx * dx + dy * dy);
        
        if (distance <= detection) {
          population(i_ind, 5) = 1;
          population(j_ind, 5) = 1;
          num_encounter[i_time] ++;
          break;
        }
      }
    }
    
    // Update individual movement pattern
    for(i_ind = 0; i_ind < num_ind; i_ind++) {
      if(population(i_ind, 5) == 1){ continue; }
      
      // speed update
      u = R::runif(0, 1) - 0.5;
      sign_u = std::copysign(1.0, u);
      abs_u = std::abs(u);
      laplace_value = acc_inter + acc_slope * population(i_ind, 3) - acc_scale * sign_u * std::log(1 - 2 * abs_u);
      population(i_ind, 3) += laplace_value;
      
      // angle update
      sign_u = std::copysign(1.0, u);
      abs_u = std::abs(u);
      population(i_ind, 2) -= turn_scale * sign_u * std::log(1 - 2 * abs_u);
      
    }
    
  }
  
  //return(num_encounter);
  return List::create(Named("encounter") = num_encounter, 
                      Named("population") = population);
  
  
}

// [[Rcpp::export]]
List one_simulation_plot(double all_L,    // field size
                             double Light_L,  // light trapped area size
                             double detection,
                             int end_time,
                             int num_ind,
                             List solo_param) {
  int i_ind, i_time, j_ind;
  double dx, dy, distance;
  double u, sign_u, abs_u, laplace_value; // for laplace random number
  
  double acc_inter  = solo_param["acc_inter"];
  double acc_slope  = solo_param["acc_slope"];
  double acc_scale  = solo_param["acc_scale"];
  double turn_scale = solo_param["turn_scale"];
  
  NumericMatrix X_time_development(num_ind, end_time + 1);
  NumericMatrix Y_time_development(num_ind, end_time + 1);
  NumericMatrix State_time_development(num_ind, end_time + 1);
  
  // Initialization
  NumericMatrix population(num_ind, 6); // Columns: X, Y, Angle, Speed, Sex, Encounter
  for(i_ind = 0; i_ind < num_ind; i_ind++) {
    population(i_ind, 0) = R::runif(0, Light_L);           // X
    population(i_ind, 1) = R::runif(0, Light_L);           // Y
    population(i_ind, 2) = R::runif(0, 2 * M_PI);          // Angle
    population(i_ind, 3) = R::runif(0, 8);                 // Speed
    population(i_ind, 4) = (i_ind < num_ind / 2) ? 0 : 1;  // Sex: 0 = female, 1 = male
    population(i_ind, 5) = 0;                              // Encounter: 0 = no encounter
    X_time_development(i_ind, 0) = population(i_ind, 0);
    Y_time_development(i_ind, 0) = population(i_ind, 1);
    State_time_development(i_ind, 0) = 0;
  }  
  

  NumericVector num_encounter(end_time);
  
  // Simulation
  for(i_time = 0; i_time < end_time; i_time++) {
    
    // Move individuals
    for(i_ind = 0; i_ind < num_ind; i_ind++) {
      if(population(i_ind, 5) == 1){ continue; }
      // Update position based on speed and angle
      population(i_ind, 0) += population(i_ind, 3) * cos(population(i_ind, 2)); // Update X
      population(i_ind, 1) += population(i_ind, 3) * sin(population(i_ind, 2)); // Update Y
      
      // Ensure individuals stay within the boundaries using wrapping
      if (population(i_ind, 0) > all_L) population(i_ind, 0) -= all_L;
      if (population(i_ind, 1) > all_L) population(i_ind, 1) -= all_L;
      if (population(i_ind, 0) < 0) population(i_ind, 0) += all_L;
      if (population(i_ind, 1) < 0) population(i_ind, 1) += all_L;
    }
    
    // Check for encounters (male and female within detection range)
    for(i_ind = 0; i_ind < num_ind / 2; i_ind++) {  // Loop over females
      if(population(i_ind, 5) == 1){ continue; }
      for (j_ind = num_ind / 2; j_ind < num_ind; j_ind++) {  // Loop over males
        if(population(j_ind, 5) == 1){ continue; }
        dx = population(i_ind, 0) - population(j_ind, 0);
        dy = population(i_ind, 1) - population(j_ind, 1);
        if (abs(dx) > all_L / 2) { dx = all_L - abs(dx); }
        if (abs(dy) > all_L / 2) { dy = all_L - abs(dy); }
        distance = sqrt(dx * dx + dy * dy);
        
        if (distance <= detection) {
          population(i_ind, 5) = 1;
          population(j_ind, 5) = 1;
          num_encounter[i_time] ++;
          break;
        }
      }
    }
    
    // Update individual movement pattern
    for(i_ind = 0; i_ind < num_ind; i_ind++) {
      if(population(i_ind, 5) == 1){ continue; }
      
      // speed update
      u = R::runif(0, 1) - 0.5;
      sign_u = std::copysign(1.0, u);
      abs_u = std::abs(u);
      laplace_value = acc_inter + acc_slope * population(i_ind, 3) - acc_scale * sign_u * std::log(1 - 2 * abs_u);
      population(i_ind, 3) += laplace_value;
      
      // angle update
      sign_u = std::copysign(1.0, u);
      abs_u = std::abs(u);
      population(i_ind, 2) -= turn_scale * sign_u * std::log(1 - 2 * abs_u);
      
    }
    
    for(i_ind = 0; i_ind < num_ind; i_ind++) {
      X_time_development(i_ind, i_time + 1) = population(i_ind, 0);
      Y_time_development(i_ind, i_time + 1) = population(i_ind, 1);
      State_time_development(i_ind, i_time + 1) = population(i_ind, 5);
    } 
  }
  
  return List::create(Named("X") = X_time_development, 
                      Named("Y") = Y_time_development, 
                      Named("State") = State_time_development);
  
}
