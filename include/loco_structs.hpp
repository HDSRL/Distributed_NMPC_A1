#ifndef LOCO_STRUCTS_H
#define LOCO_STRUCTS_H

#include <Eigen/Dense>
#include <casadi/casadi.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <filesystem>

#define SIM_TIME 30.0
#define HORIZON 15
#define GOAL_X 10.0
#define GOAL_Y 0.0
#define NUMBER_OF_AGENTS 2
#define NUMBER_OF_OBS 5
// #define GOAL_X 10.0
// #define GOAL_Y 0.0
// #define NUMBER_OF_AGENTS 1

using namespace casadi;
namespace fs = std::filesystem;

SX blkdiag(const std::vector<SX>& matrices);

struct Params {
    size_t mpc_iter = 0; // Iteration counter
    double sim_length = 10.0; // Default simulation time
    double current_time = 0.0; // Initialize to zero
    double dt; // Sample time
    double Fs; // Control Frequency, initialized in constructor
    double mu = 0.6; // Friction Coefficient
    size_t L; // Simulation length, initialized in constructor
    size_t N = HORIZON; // Horizon length
    double height = 0.28; // Stand height

    size_t n_states = 12; // Number of states
    size_t n_inputs = 12; // Number of inputs
    size_t n_controls = 12;
    size_t gait = 8; // Trot

    double G = 9.81; // Gravity
    double mass = 12.4530; // Mass

    DM J = DM::zeros(3,3); // Assuming a 3x3 inertia matrix

    // Nominal hip position for most gaits
    DM pf = DM::zeros(3, 4); 
    DM p_hip = DM::zeros(3, 4);

    Eigen::MatrixXd Pr_refined_; // To be initialized as needed
    Eigen::MatrixXd Prd_refined_; // To be initialized as needed
    size_t agentIndex = 1;

    DM contact = DM::vertcat({1, 1, 1, 1});
    DM desVel = DM::vertcat({1, 0});

    SX Q; // State cost matrix
    SX R; // Input cost matrix

    size_t stepcnt = 0;

    // Constructor to initialize members that depend on dt
    explicit Params() {
        
        SX Qp = diag(SX::vertcat({8e8, 8e8, 8e8}));
        SX Qv = diag(SX::vertcat({1e3, 1e3, 1e3}));
        SX Qr = diag(SX::vertcat({8e4, 8e4, 3e3}));
        SX Qw = 1.0 * SX::eye(3);

        // Combine into a block diagonal matrix Q
        Q = blkdiag({Qp, Qv, Qr, Qw}); // Corrected by adding braces to create a vector

        // Input cost component and R matrix
        SX R_tmp = 1e-2 * SX::eye(3);
        R = blkdiag({R_tmp, R_tmp, R_tmp, R_tmp}); // Corrected by adding braces to create a vector    
        
        // Initialize any other members that require dynamic calculation here
    }
};

struct Args {
    DM p, x0, lbx, ubx, lbg, ubg; 
};

#endif // LOCO_STRUCTS_H