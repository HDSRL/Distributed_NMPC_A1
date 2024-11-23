#ifndef NMPC_C
#define NMPC_C

#include "loco_structs.hpp"
#include <chrono>
#include <yaml-cpp/yaml.h>

using namespace casadi;
namespace fs = std::filesystem;

// #define NUMBER_OF_AGENTS 1
// #define NUMBER_OF_OBS 9

class NMPC {

private:
    
    int cbf_active = 0;
    Args args;     // Arguments for the solver
    Function f;
    double curr_time;
    Function f_dyn_;
    Function solver_;
    std::string basePath;
    bool isDomainChanged_ = false;

    // NMPC variables
    DM qk_; // Current state (qk) in the other code
    DM qkAll_; // All states for N horizons (X0) in the other code
    int domain_ = 0;
    int control_tick_ = 0;

    DM q_des_; // Desired state (xd) in the other code
    DM q_des_log_; // Desired state (xd) in the other code
    DM uk_; // Current control input - n_inputs x 1 
    DM ukAll_; // All control inputs for N horizons (u0) in the other code
    DM x_sol_; // Current state (xk) in the other code
    DM x_sol_N_; // solution for all states in shape n_states x N
    DM x_sol_N_w_IC; // solution for all states in shape n_states x N

    DM bez_alpha_; // Bezier coefficients for COM trajectory
    DM bez_GRF_; // Bezier coefficients for GRF trajectory
    
    // Logger variables
    DM uout;
    DM xout;
    DM feet;
    DM traj;
    DM tout;
    DM cbfout;
    DM NMPC_solve_time;

    // From low-level controller
    double q[18] = {0};
    double dq[18] = {0};
    Eigen::Vector4d state_other_;  // OSDCP state data := {x, y, dx, dy}
    int contactInd[4] = {1};
    Eigen::Matrix<double, 3, 4> toePos_;

    // Solver inputs
    DM g_ub_, g_lb_;

    // Eigen vars
    Eigen::MatrixXd x_sol_eigen_; // mpc_state_e_x_eventbased_ in the old code - this is only state solution for N domains
    Eigen::MatrixXd GRF_sol_eigen_;
    Eigen::MatrixXd GRF_sol_other_;
    Eigen::MatrixXd bez_alpha_eigen_; // alpha_COM_traj_e_ in the old code
    Eigen::MatrixXd bez_GRF_eigen_;
    Eigen::Matrix<double, 2*NUMBER_OF_AGENTS, 1> Pstart_;
    Eigen::Matrix<double, 2,NUMBER_OF_OBS> Pobs;
    Eigen::Matrix<double, 2,NUMBER_OF_OBS> Pobs_real;

    size_t agent_id_;
    Eigen::Vector2d agent_Initial_;
    
    Eigen::MatrixXd Pr_refined_;
    Eigen::MatrixXd Prd_refined_;

    Eigen::Vector2d obs_near_;
    YAML::Node config;

    std::ofstream xoutFile_;
    std::ofstream trajFile_;
    std::ofstream uoutFile_;
    std::ofstream feetFile_;
    std::ofstream toutFile_;
    std::ofstream cbfFile_;
    std::ofstream nmpcSolveTimeFile_;

public:
    Params params;

    NMPC()
    {
        // Initialize the parameters
        params.mpc_iter = 0;
        params.current_time = 0.0;
        params.dt = 0.01;
        params.Fs = 1/params.dt;
        params.sim_length = SIM_TIME;
        // params.L = 1000; // not sure if this is needed
        params.N = HORIZON;   // Horizon length
        params.height = 0.28;   // stand height
        params.n_states = 12;   // states
        params.n_inputs = 12;   // control inputs
        params.n_controls = 12; // control inputs
        params.gait = 8;        // Trot gait
        params.G = 9.81;        // gravity
        params.mass = 12.4530;  // mass
        params.J = DM::zeros(3,3);  // inertia matrix
        params.pf = DM::zeros(3, 4);    // nominal foot position matrix
        params.p_hip = DM::zeros(3, 4); // nominal hip position matrix
        params.Pr_refined_ = Eigen::MatrixXd::Zero(4, 4);
        params.Prd_refined_ = Eigen::MatrixXd::Zero(4, 4);
        params.agentIndex = 1;  // agent index

        // Initialize the pf matrix within the constructor body
        params.pf(0,0) = 0.183;  params.pf(1,0) = -0.1321; params.pf(2,0) = 0.01;
        params.pf(0,1) = 0.183;  params.pf(1,1) = 0.1321;  params.pf(2,1) = 0.01;
        params.pf(0,2) = -0.183; params.pf(1,2) = -0.1321; params.pf(2,2) = 0.01;
        params.pf(0,3) = -0.183; params.pf(1,3) = 0.1321;  params.pf(2,3) = 0.01;

        params.p_hip(0,0) = 0.183;  params.p_hip(1,0) = -0.1321; params.p_hip(2,0) = 0.01;
        params.p_hip(0,1) = 0.183;  params.p_hip(1,1) = 0.1321;  params.p_hip(2,1) = 0.01;
        params.p_hip(0,2) = -0.183; params.p_hip(1,2) = -0.1321; params.p_hip(2,2) = 0.01;
        params.p_hip(0,3) = -0.183; params.p_hip(1,3) = 0.1321;  params.p_hip(2,3) = 0.01;

        
        params.J(0,0) = 0.01683993;  params.J(0,1) = 8.3902e-5;   params.J(0,2) = 0.000597679;
        params.J(1,0) = 8.3902e-5;   params.J(1,1) = 0.056579028; params.J(1,2) = 2.5134e-5;
        params.J(2,0) = 0.000597679; params.J(2,1) = 2.5134e-5;   params.J(2,2) = 0.064713601;

        basePath = "tmp/";

        config = YAML::LoadFile("../src/config.yaml");
    };
    ~NMPC();

    casadi::DM bezier(const DM& afra, double s);
    casadi::DM bezierd(const DM& alpha, double s);

    casadi::SX Rotz(const casadi::SX& t);
    casadi::SX skewSym(const casadi::SX& x);
    void planner();
    void force_inequality(const casadi::SX& U, casadi::SX& g_force);
    void writeMatrixToFile(const casadi::SX& matrix, const std::string& filename);
    void writeMatrixToFile(const casadi::DM& matrix, const std::string& filename);
    void writeMatrixToFile(const std::string& filename);
    void appendMatrixToFile(const casadi::DM& matrix, std::ofstream& file);
    void equalityconstraints(const casadi::SX& X, const casadi::SX& U, const casadi::SX& P, casadi::SX& obj, casadi::SX& g_eq, const Function& f);
    void set_constraint_bounds();
    void set_statebounds();
    void force_inequality_bounds();
    void set_collision_bounds();
    void shift();
    int binomial(int n, int k);
    int binomialCoeff(int n, int k);
    void fitComTrajectory();
    void fitGRFTrajectory();
    void casadiToEigen(const casadi::DM& mat1, Eigen::MatrixXd& mat2);
    void eigenToCasadi(const Eigen::MatrixXd& vec, casadi::DM& mat);
    void eigenToCasadi(const Eigen::Vector4d& vec, casadi::DM& mat);
    void updateState(double* qin, double* dqin, int* ind, Eigen::Matrix<double, 3, 4> toePos, Eigen::Vector4d state_vec, Eigen::VectorXd GRF_sol_other, const int domain, const int control_tick);
    void generateReferenceTrajectory();
    void find_nearest_obs();

    void setPstart(const Eigen::Matrix<double, 2*NUMBER_OF_AGENTS, 1>& Pstart);
    void setPobs(Eigen::Matrix<double, 2,NUMBER_OF_OBS>& Pobs);
    void setPobs_real(Eigen::Matrix<double, 2,NUMBER_OF_OBS>& Pobs);
    void setAgentID(size_t agent_id);
    int findNearestIndex();

    void code_gen();
    void run_NMPC();
    void init();
    void solve_nlp();
    void logData();
    inline std::vector<double> get_com_pos() {
        casadi::DM first_three = qk_(casadi::Slice(0, 3));  // Extract first three elements
        return std::vector<double>(first_three->begin(), first_three->end());
    }

    // Modify this function to return std::vector<double> for uk_
    inline std::vector<double> get_GRF() {
        return std::vector<double>(uk_->begin(), uk_->end());  // Convert the whole uk_ to std::vector
    }

    inline std::vector<double> get_contacts() {
        return std::vector<double>(params.contact->begin(), params.contact->end());  // Convert the whole uk_ to std::vector
    }

    inline std::vector<double> get_Pf() {
        casadi::DM Pf = casadi::DM::vertcat({
            params.pf(Slice(), 0), // First column
            params.pf(Slice(), 1), // Second column
            params.pf(Slice(), 2), // Third column
            params.pf(Slice(), 3)  // Fourth column
        });
        return std::vector<double>(Pf->begin(), Pf->end());  // Convert the whole uk_ to std::vector
    }

    inline Eigen::MatrixXd get_alphaCOM() {return bez_alpha_eigen_;}
    inline Eigen::MatrixXd get_alphaGRF() {return bez_GRF_eigen_;}
    inline Eigen::MatrixXd get_MPCsol() {return x_sol_eigen_;}
    inline Eigen::MatrixXd get_GRFsol() { return GRF_sol_eigen_; }
    inline int getCBF_status() { return cbf_active; }
    Eigen::Vector4d get_lastState();
};
#endif // NMPC_H