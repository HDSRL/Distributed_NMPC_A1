#ifndef NMPC_P
#define NMPC_P

#include <Eigen/Dense>
#include <casadi/casadi.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace casadi;

SX blkdiag(const std::vector<SX>& matrices);
void writeMatrixToFile(const casadi::DM& matrix, const std::string& filename);

struct Params {
    size_t mpc_iter = 0; // Iteration counter
    double current_time = 0.0; // Initialize to zero
    double dt; // Sample time
    double Fs; // Control Frequency, initialized in constructor
    double mu = 0.6; // Friction Coefficient
    size_t L; // Simulation length, initialized in constructor
    size_t N; // Horizon length
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
    DM desVel = DM::vertcat({1.0, 0});

    SX Q; // State cost matrix
    SX Q_terminal; // Terminal cost matrix
    SX R; // Input cost matrix

    size_t stepcnt = 0;

    // Constructor to initialize members that depend on dt
    explicit Params(double sampleTime) : dt(sampleTime), Fs(1.0 / dt), L(static_cast<int>(20 * Fs)), N(15) {
        // Initialize the pf matrix within the constructor body
        pf(0,0) = 0.183;  pf(1,0) = -0.1321; pf(2,0) = 0.01;
        pf(0,1) = 0.183;  pf(1,1) = 0.1321;  pf(2,1) = 0.01;
        pf(0,2) = -0.183; pf(1,2) = -0.1321; pf(2,2) = 0.01;
        pf(0,3) = -0.183; pf(1,3) = 0.1321;  pf(2,3) = 0.01;

        p_hip(0,0) = 0.183;  p_hip(1,0) = -0.1321;  p_hip(2,0) = 0.01;
        p_hip(0,1) = 0.183;  p_hip(1,1) =  0.1321;  p_hip(2,1) = 0.01;
        p_hip(0,2) = -0.183; p_hip(1,2) = -0.1321;  p_hip(2,2) = 0.01;
        p_hip(0,3) = -0.183; p_hip(1,3) =  0.1321;  p_hip(2,3) = 0.01;

        
        J(0,0) = 0.01683993;  J(0,1) = 8.3902e-5;   J(0,2) = 0.000597679;
        J(1,0) = 8.3902e-5;   J(1,1) = 0.056579028; J(1,2) = 2.5134e-5;
        J(2,0) = 0.000597679; J(2,1) = 2.5134e-5;   J(2,2) = 0.064713601;

        SX Qp = diag(SX::vertcat({1e8, 1e9, 1e2}));
        SX Qv = diag(SX::vertcat({1e8, 1e9, 1e2}));
        SX Qr = diag(SX::vertcat({1e8, 1e9, 1e-1}));
        SX Qw = diag(SX::vertcat({1e8, 1e9, 1e-1}));//1.0 * SX::eye(3);

        SX Qp_terminal = diag(SX::vertcat({7e6, 1e6, 1e4}));
        SX Qv_terminal = diag(SX::vertcat({7e6, 1e6, 1e4}));
        SX Qr_terminal = diag(SX::vertcat({1e3, 1e3, 1e-1}));
        SX Qw_terminal = diag(SX::vertcat({1e3, 1e3, 1e-1}));

        // Combine into a block diagonal matrix Q
        Q = blkdiag({Qp, Qv, Qr, Qw}); // Corrected by adding braces to create a vector
        Q_terminal = blkdiag({Qp_terminal, Qv_terminal, Qr_terminal, Qw_terminal});

        // Input cost component and R matrix
        SX R_tmp = 1e-2 * SX::eye(3);
        R = blkdiag({R_tmp, R_tmp, R_tmp, R_tmp}); // Corrected by adding braces to create a vector    
        
        // Initialize any other members that require dynamic calculation here
    }
};

SX blkdiag(const std::vector<SX>& matrices) {
    if (matrices.empty()) return SX();

    // Start with the first matrix
    SX result = matrices[0];
    for (size_t i = 1; i < matrices.size(); ++i) {
        // Pad the existing result and the next matrix with zeros to make them block diagonal
        SX topRight = SX::zeros(result.size1(), matrices[i].size2());
        SX bottomLeft = SX::zeros(matrices[i].size1(), result.size2());

        // Concatenate to form the new result
        result = vertcat(horzcat(result, topRight), 
                         horzcat(bottomLeft, matrices[i]));
    }

    return result;
}

struct Args {
    DM p, x0, lbx, ubx, lbg, ubg; 
};

DM bezier(const DM& afra, double s) {
    int n = afra.size1(); // Number of rows
    int m = afra.size2(); // Number of columns
    DM value = DM::zeros(n);
    int M = m - 1;

    std::vector<int> k; // Binomial coefficients
    switch (M) {
        case 3: k = {1, 3, 3, 1}; break;
        case 4: k = {1, 4, 6, 4, 1}; break;
        case 5: k = {1, 5, 10, 10, 5, 1}; break;
        case 6: k = {1, 6, 15, 20, 15, 6, 1}; break;
        case 7: k = {1, 7, 21, 35, 35, 21, 7, 1}; break;
        case 8: k = {1, 8, 28, 56, 70, 56, 28, 8, 1}; break;
        case 9: k = {1, 9, 36, 84, 126, 126, 84, 36, 9, 1}; break;
        case 10: k = {1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1}; break;
        case 20: k = {1, 20, 190, 1140, 4845, 15504, 38760, 77520, 125970, 167960, 184756, 167960, 125970, 77520, 38760, 15504, 4845, 1140, 190, 20, 1}; break;
        default: std::cerr << "Degree " << M << " not supported." << std::endl; return value; // Optionally handle other degrees dynamically
    }

    std::vector<double> x(M + 1, 1.0), y(M + 1, 1.0);
    for (int i = 1; i <= M; ++i) {
        x[i] = s * x[i - 1];
        y[i] = (1 - s) * y[i - 1];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= M; ++j) {
            value(i) += afra(i, j) * k[j] * x[j] * y[M - j];
        }
    }

    return value;
}

DM bezierd(const DM& alpha, double s) {
    int n = alpha.size1(); // Number of rows
    int m = alpha.size2(); // Number of columns
    DM value = DM::zeros(n);
    int M = m - 1; // Degree of the Bezier curve

    std::vector<int> k; // Binomial coefficients for the derivative
    switch (M) {
        case 3: k = {3, 6, 3}; break;
        case 4: k = {4, 12, 12, 4}; break;
        case 5: k = {5, 20, 30, 20, 5}; break;
        case 6: k = {6, 30, 60, 60, 30, 6}; break;
        case 7: k = {7, 42, 105, 140, 105, 42, 7}; break;
        case 8: k = {8, 56, 168, 280, 280, 168, 56, 8}; break;
        case 9: k = {9, 72, 252, 504, 630, 504, 252, 72, 9}; break;
        case 10: k = {10, 45, 120, 210, 252, 210, 120, 45, 10, 1}; break;
        case 20: k = {20, 380, 3420, 19380, 77520, 232560, 542640, 1007760, 1511640, 1847560, 1847560, 1511640, 1007760, 542640, 232560, 77520, 19380, 3420, 380, 20}; break;
        default: std::cerr << "Degree " << M << " not supported for derivative calculation." << std::endl; return value;
    }

    std::vector<double> x(M, 1.0), y(M, 1.0); // Adjust for derivative size
    for (int i = 1; i < M; ++i) {
        x[i] = s * x[i - 1];
        y[i] = (1 - s) * y[i - 1];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < M; ++j) {
            value(i) += (alpha(i, j + 1) - alpha(i, j)) * k[j] * x[j] * y[M - 1 - j];
        }
    }

    return value;
}

SX Rotz(const SX& t) {
    // Ensure t is a vector of size 3
    if (t.size1() != 3 || t.size2() != 1) {
        throw std::runtime_error("Input vector t must be of size 3.");
    }

    casadi::SX Rz = casadi::SX::zeros(3, 3);

    // Note: In CasADi, indexing is zero-based, similar to C++ and Eigen but different from MATLAB.
    Rz(0, 0) = cos(t(1)) * cos(t(2));
    Rz(0, 1) = -sin(t(2));
    Rz(0, 2) = 0;
    Rz(1, 0) = cos(t(1)) * sin(t(2));
    Rz(1, 1) = cos(t(2));
    Rz(1, 2) = 0;
    Rz(2, 0) = -sin(t(1));
    Rz(2, 1) = 0;
    Rz(2, 2) = 1;

    return Rz;
}

casadi::SX skewSym(const casadi::SX& x) {
    using namespace casadi;

    // Determine the size based on the input x
    int cols = x.size2();
    SX mat = SX::zeros(3, 3 * cols);

    // Loop through each column of x to create skew-symmetric matrices
    for (int n = 0; n < cols; ++n) {
        SX slice = x(Slice(), n);  // Get the nth column of x
        SX skew = SX::zeros(3, 3);
        // Fill the skew-symmetric matrix
        skew(0, 1) = -slice(2);
        skew(0, 2) = slice(1);
        skew(1, 0) = slice(2);
        skew(1, 2) = -slice(0);
        skew(2, 0) = -slice(1);
        skew(2, 1) = slice(0);

        // Place the skew matrix into the correct position in mat
        mat(Slice(0, 3), Slice(3*n, 3*(n+1))) = skew;
    }

    return mat;
}

void planner(double t, const DM& q, DM& q_des, int mode, Params& params) {
    const int n_states = 12; // Example: Number of states
    // const int N = params.N;
    DM xd = DM::zeros(n_states * params.N);

    if (t <= 1.0) { // Stand up duration
        mode = 1;
    }

    DM traj = q_des(casadi::Slice(0, n_states), 0); // Second parameter 0 is to indicate column index for a vector


    switch (mode) {
        case 1: {
            double t0 = 0, tf = 1.0; // Start and end times for mode 1
            double s = std::min((t - t0) / (tf - t0), 1.0); // Ensure s is within [0, 1]
            
            // Extract relevant parameters for easy access
            double xf = 0.0; // Initial x position
            double yf = 0.0; // Initial y position
            double zf = params.height; // Assuming the height is constant or predefined
            double z0 = 0.08;

            // Using DM::horzcat to create each row of the alpha matrix
            casadi::DM row1 = casadi::DM::horzcat({DM(0), DM(0), DM(0), DM(xf / 4), DM(3 * xf / 4), DM(xf), DM(xf), DM(xf)});
            casadi::DM row2 = casadi::DM::horzcat({DM(0), DM(0), DM(0), DM(yf / 4), DM(3 * yf / 4), DM(yf), DM(yf), DM(yf)});
            casadi::DM row3 = casadi::DM::horzcat({DM(z0), DM(z0), DM(z0), DM(zf / 4), DM(3 * zf / 4), DM(zf), DM(zf), DM(zf)});

            casadi::DM alpha = casadi::DM::vertcat({row1, row2, row3});

            // Update trajectory position
            casadi::DM position_update = bezier(alpha, s);

            for (int i = 0; i < 3; ++i) {
                traj(i) = position_update(i);
            }

            // Update trajectory velocity
            casadi::DM velocity_update = bezierd(alpha, s);

            for (int i = 0; i < 3; ++i) {
                traj(i + 3) = velocity_update(i);
            } 
            break;
        }
        case 8: {
            double domLen = 200.0; // Domain length in milliseconds
            double domLenSec = domLen / 1000.0; // Domain length in seconds
            DM Rz = Rotz(vertcat(q(6), q(7), q(8))); // Assuming q contains Euler angles in the last 3 elements
            // size_t currentStep = std::floor(t / params.dt) + 1;

            if (t >= (params.stepcnt * domLenSec + 1)) {
                // Update contact sequence
                if (static_cast<double>(params.contact(0)) == 1.0) {
                    params.contact = {0, 1, 1, 0};
                } else {
                    params.contact = {1, 0, 0, 1};
                }

                // Increase velocity in body frame
                traj(3) += 0.1; // Increase command velocity in body X
                traj(4) += 0.1; // Increase command velocity in body Y

                // Saturate the command velocities
                traj(3) = fmin(traj(3), params.desVel(0));
                traj(4) = fmin(traj(4), params.desVel(1));

                // Calculate step length in body frame and update step foot position
                DM stepLen = vertcat(traj(3) * domLenSec / 2, traj(4) * domLenSec / 2, 0.0);
                DM hipTmp = params.p_hip + stepLen;
                hipTmp = mtimes(Rz, hipTmp); // Rotate to the world frame
                params.pf = hipTmp + vertcat(q(Slice(0, 2)), 0.0);

                // Set velocity in world frame
                DM cmdWrld = mtimes(Rz, traj(Slice(3, 6))); // Convert to world frame

                // std::string basePath = "tmp/pathplanner/";
                // writeMatrixToFile(stepLen, basePath + "stepLen.txt");
                // writeMatrixToFile(params.p_hip, basePath + "params_p_hip.txt");
                // writeMatrixToFile(hipTmp, basePath + "hipTmp.txt");
                // writeMatrixToFile(q(Slice(0, 3)), basePath + "q_COM.txt");
                // writeMatrixToFile(cmdWrld, basePath + "cmdWrld.txt");
                // writeMatrixToFile(params.pf, basePath + "params_pf.txt");



                traj(Slice(3, 6)) = cmdWrld;

                // Increment step counter
                params.stepcnt++;
            }

            float heading = 0;

            // Update position and height based on velocity
            traj(0) = q(0) + params.dt * traj(3);
            traj(1) = q(1) + params.dt * traj(4);
            traj(2) = params.height; // Assuming height is constant or predefined
            traj(5) = 0; // Assuming zero vertical velocity

            for (int n = 1; n <= params.N; ++n) {
                // Compute the predicted state for each step in the horizon
                DM state = DM::zeros(n_states, 1);
                
                // Updating the state based on your MATLAB code logic
                state(0) = q(0) + n * params.dt * traj(3); // pos X
                state(1) = q(1) + n * params.dt * traj(4); // pos Y
                state(2) = params.height;                  // pos Z (constant height assumption)
                state(3) = traj(3);                        // vel X (constant velocity assumption)
                state(4) = traj(4);                        // vel Y
                state(5) = 0.0;                            // vel Z (assuming zero vertical velocity)
                // Assuming Euler angles (heading) and angular velocity are zero
                state(6) = 0.0; // Roll or phi
                state(7) = 0.0; // Pitch or theta
                state(8) = heading; // Yaw or psi
                state(9) = 0.0;  // Angular velocities (assuming zeros)
                state(10) = 0.0;
                state(11) = 0.0;

                // Place the computed state into the appropriate segment of xd
                //xd.segment((n-1)*n_states, n_states) = state;
                xd(Slice((n-1)*n_states, n*n_states)) = state;
            }
            break;
        }
        // Handle other cases as necessary
        default:
            std::cerr << "Mode '" << mode << "' is not defined in the motion planner.\n";
    }

    if (mode != 8) {
        xd = repmat(traj, params.N, 1);
    }

    q_des = xd; // Only xd needs to be returned
}

void force_inequality(const casadi::SX& U, const Params& params, casadi::SX& g_force, casadi::SX& g_ub, casadi::SX& g_lb) {
    // Use SX for symbolic constraint expressions
    casadi::SX g_force_sx = casadi::SX::vertcat(std::vector<casadi::SX>{});

    // g_ub and g_lb as symbolic vectors with appropriate bounds
    casadi::SX g_ub_sx = casadi::SX::vertcat(std::vector<casadi::SX>{});
    casadi::SX g_lb_sx = casadi::SX::vertcat(std::vector<casadi::SX>{});


    for (size_t k = 0; k < params.N; ++k) {
        for (size_t leg = 0; leg < 4; ++leg) {
            size_t base_idx = leg * 3;
            // Directly construct and append symbolic expressions
            g_force_sx = casadi::SX::vertcat({g_force_sx,
                U(base_idx, k) - params.mu * U(base_idx + 2, k),
                U(base_idx + 1, k) - params.mu * U(base_idx + 2, k),
                U(base_idx, k) + params.mu * U(base_idx + 2, k),
                U(base_idx + 1, k) + params.mu * U(base_idx + 2, k)
            });
            g_ub_sx = casadi::SX::vertcat({g_ub_sx, 0.0, 0.0, casadi::inf, casadi::inf}); 
            g_lb_sx = casadi::SX::vertcat({g_lb_sx, -casadi::inf, -casadi::inf, 0, 0});
        }
    }


    // Convert symbolic SX expressions to DM for g_ub and g_lb if necessary
    // Note: This might not be needed if keeping them symbolic is preferable
    g_force = g_force_sx; // g_force is already SX and can stay symbolic
    g_ub = g_ub_sx;
    g_lb = g_lb_sx;
}

void force_inequality_bounds(const Params& params, casadi::DM& g_ub, casadi::DM& g_lb) {

    // g_ub and g_lb as symbolic vectors with appropriate bounds
    casadi::DM g_ub_sx;
    casadi::DM g_lb_sx;

    for (size_t k = 0; k < params.N; ++k) {
        for (size_t leg = 0; leg < 4; ++leg) {
            size_t base_idx = leg * 3;
            // Directly construct and append symbolic expressions
            g_ub_sx = casadi::SX::vertcat({g_ub_sx, 0.0, 0.0, casadi::inf, casadi::inf}); 
            g_lb_sx = casadi::SX::vertcat({g_lb_sx, -casadi::inf, -casadi::inf, 0, 0});
        }
    }
    g_ub = g_ub_sx;
    g_lb = g_lb_sx;
}

// Function to write casadi::SX matrix to file
void writeMatrixToFile(const casadi::SX& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Print matrix dimensions
        // file << "Matrix size: " << matrix.size1() << "x" << matrix.size2() << std::endl;
        // Print matrix content preserving original shape
        for (int i = 0; i < matrix.size1(); ++i) {
            for (int j = 0; j < matrix.size2(); ++j) {
                file << matrix(i, j);
                if (j < matrix.size2() - 1) file << ", ";
            }
            file << std::endl;
        }
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

// Function to write casadi::DM matrix to file
void writeMatrixToFile(const casadi::DM& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        // Print matrix dimensions
        // file << "Matrix size: " << matrix.size1() << "x" << matrix.size2() << std::endl;
        // Print matrix content preserving original shape
        for (int i = 0; i < matrix.size1(); ++i) {
            for (int j = 0; j < matrix.size2(); ++j) {
                file << matrix(i, j);
                if (j < matrix.size2() - 1) file << ", ";
            }
            file << std::endl;
        }
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

// Overload to handle multiple DM matrices and one size_t variable
void writeMatrixToFile(const Params& params, const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::trunc);  // Open for output and truncate to clear existing content
    if (file.is_open()) {
        // Helper lambda to write a single matrix
        auto writeMatrix = [&](const casadi::DM& matrix) {
            file << "Matrix size: " << matrix.size1() << "x" << matrix.size2() << "\n" << std::endl;
            for (int i = 0; i < matrix.size1(); ++i) {
                for (int j = 0; j < matrix.size2(); ++j) {
                    file << matrix(i, j);
                    if (j < matrix.size2() - 1) file << ", ";
                }
                file << std::endl;
            }
        };

        // Write each matrix with a header
        file << "params.pf:" << std::endl;
        writeMatrix(params.pf);
        file << "\nparams.contact:" << std::endl;
        writeMatrix(params.contact);
        file << "\nparams.desVel:" << std::endl;
        writeMatrix(params.desVel);
        
        // Additionally write the size_t variable
        file << "\nparams.stepcnt: " << params.stepcnt << std::endl;
        file << "\n\nparams.current_time: " << params.current_time << std::endl;
        file << "\n\nparams.mpc_iter: " << params.mpc_iter << std::endl;

        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
}

void equalityconstraints(const SX& X, const SX& U, const SX& P, const Params& params, SX& obj, SX& g_eq, SX& slack_var, const Function& f) {

    std::string basePath = "tmp/";
    int n_states = params.n_states;
    int N = params.N;
    double dt = params.dt;

    // Initial condition constraints
    SX st = X(Slice(0, n_states), 0);

    int end_index = P.numel();
    SX pf = P(Slice(end_index - 12, end_index), 0); // Slice the last 12 elements

    g_eq = st - P(Slice(0, n_states));


    double sl_gain = 5e9;
    // TODO - load it in real-time from yaml file 
    obj = sl_gain*mtimes(slack_var.T(), slack_var);

    SX st_x_avg = 0;
    SX st_y_avg = 0;

    for (int k = 0; k < params.N; ++k) {
        SX st = X(Slice(0, n_states), k);
        SX con = U(Slice(0, n_states), k);

        SX des_st = P(Slice(k * n_states, (k + 1) * n_states)); // Desired state

        // writeMatrixToFile(st, basePath + "st.txt");
        // writeMatrixToFile(des_st, basePath + "des_st.txt");
        // writeMatrixToFile(params.Q, basePath + "params_Q.txt");
        // writeMatrixToFile(con, basePath + "con.txt");
        // writeMatrixToFile(params.R, basePath + "params_R.txt");
        // Calculate objective
        SX st_err = st - des_st;
        // writeMatrixToFile(st_err, basePath + "st_err.txt");
        if (k < N)
        {
            SX stage_st_cost = mtimes(mtimes(st_err.T(), params.Q), st_err);
            SX stage_con_cost = mtimes(mtimes(con.T(), params.R), con);
            obj = obj + stage_st_cost + stage_con_cost;
        }
        else
        {
            SX terminal_st_cost = mtimes(mtimes(st_err.T(), params.Q_terminal), st_err);
            SX terminal_con_cost = mtimes(mtimes(con.T(), params.R), con);
            obj = obj + terminal_st_cost + terminal_con_cost;
        }        

        // writeMatrixToFile(stage_cost, basePath + "stage_cost.txt");
        // writeMatrixToFile(terminal_cost, basePath + "terminal_cost.txt");


        // Calculate next state prediction
        SX st_next = X(Slice(0, n_states), k + 1);

        std::map<std::string, SX> arg;
        arg["st"] = st;  // Example random inputs
        arg["con"] = con;
        arg["pf"] = pf;

        std::map<std::string, SX> result = f(arg);
        SX f_value = result.at("rhs");

        // writeMatrixToFile(st_next, basePath + "st_next.txt");
        // writeMatrixToFile(f_value, basePath + "f_value.txt");

        // Compute constraints
        SX st_next_euler = st + dt * f_value;
        g_eq = SX::vertcat({g_eq, st_next - st_next_euler});
        
        st_x_avg = st_x_avg + st(0);
        st_y_avg = st_y_avg + st(1);

        // writeMatrixToFile(st_next_euler, basePath + "st_next_euler.txt");
        // writeMatrixToFile(g_eq, basePath + "g_eq.txt");
    }

    // Add agreement protocol: LJ potential
    SX other_agent_st = DM::zeros(2, 1);
    int len = P.numel(); 
    other_agent_st(0) = P(len - 2 - 2 -16 -12);
    other_agent_st(1) = P(len - 1 - 2 -16 -12); // -12 for footholds which are stacked at the end

    SX dist = sqrt(pow(st_x_avg/15 - other_agent_st(0), 2) + pow(st_y_avg/15 - other_agent_st(1), 2)); 

    // obj = obj + 0.2e8*4*0.5*(pow(0.8/dist, 12) - pow(0.8/dist, 6)); 
    
}

void set_constraint_bounds(Args& args, const DM& g_lb, const DM& g_ub, int n_states, int N) {
    // Calculate the start index for setting lower and upper bounds

    int startIndex = n_states * (N + 1); // Adjust for zero-based indexing

    // Resize the lbg and ubg if necessary
    if (args.lbg.is_empty() || args.lbg.size1() < startIndex + g_lb.size1()) {
        args.lbg = DM::zeros(startIndex + g_lb.size1(), 1); // Adjust total size if needed
    }
    if (args.ubg.is_empty() || args.ubg.size1() < startIndex + g_ub.size1()) {
        args.ubg = DM::zeros(startIndex + g_ub.size1(), 1); // Adjust total size if needed
    }

    // Set values for lower bounds
    for (int i = 0; i < g_lb.size1(); ++i) {
        args.lbg(startIndex + i) = g_lb(i);
    }

    // Set values for upper bounds
    for (int i = 0; i < g_ub.size1(); ++i) {
        args.ubg(startIndex + i) = g_ub(i);
    }
}

void set_statebounds(int mpciter, const Params& params, Args& args) {
    int n_states = params.n_states;  // Number of states
    int n_inputs = params.n_inputs;  // Number of inputs
    int N = params.N;
    double zmax = 150;

    // Initialize LB and UB
    args.lbx = DM::zeros(n_states*(N+1) + n_inputs*N + 1, 1); // +1 for slack_var
    args.ubx = DM::zeros(n_states*(N+1) + n_inputs*N + 1, 1); // +1 for slack_var
    args.lbx(Slice()) = -DM::inf();
    args.ubx(Slice()) = DM::inf();

    for (int k = 0; k < N; ++k) {
        int contact_index = (mpciter + k) % 40; // Adjusting MATLAB 1-based to C++ 0-based index
        for (int leg = 0; leg < 4; ++leg) {
            int index = n_states*(N+1) + n_inputs*k + 3*leg + 2; // Adjusted for 0-based indexing
            args.lbx(index) = 0;
            args.ubx(index) = params.contact(leg) * zmax; // Assuming contact is std::vector<double>
        }
    }
}

void shift(double& t0, DM& qk, DM& u0, double dt, const DM& u, const Function& f) {
    // Prepare input for the function f
    std::map<std::string, DM> arg;
    arg["st"] = qk;  // Current state
    arg["con"] = u(0, Slice()).T();  // First control input, transposed to match MATLAB

    // Call the function f
    std::map<std::string, DM> result = f(arg);
    DM f_value = result["rhs"];  // Assuming the output is named 'rhs'

    // Update state
    qk = qk + DM(dt) * f_value;

    // Increment time
    t0 += dt;

    // Shift control inputs
    DM u_new = vertcat(u(Slice(1, u.size1()), Slice()),  // Skip the first control set
                       u(Slice(u.size1()-1, u.size1()), Slice()));  // Duplicate the last control set
    u0 = u_new;
}

void collision_constraints(const casadi::SX& X, const casadi::SX& P, const Params& params, casadi::SX& g_obs, Args& args) {
    using namespace casadi;

    casadi::SX o_lb, o_ub;
    int n_states = params.n_states; // Extract parameters
    // Extract parameters
    int N = params.N; // Make sure to convert the type correctly
    double d_safe = 0.5;


    // SX obs = P(Slice(N * n_states, (N + 1) * n_states)); // Obstacle position

    // TODO - find nearest obstacle in real-time
    // SX obs = SX::vertcat({4, -0.2}); // Obstacle position
    int end_idx_P = P.numel();
    SX obs = vertcat(P(end_idx_P-2 -16 -12), P(end_idx_P-1 -16 -12)); // Obstacle position

    // Initialize constraint array
    g_obs = SX::zeros(params.N, 1); // Use vertcat for a column vector

    // Loop over each time step
    for (int k = 0; k < N; ++k) {
        // Extract x and y positions at step k
        SX x_pos = X(0, k+1); // Adjust index for C++ (0-based)
        SX y_pos = X(1, k+1);

        // Calculate Euclidean distance from the obstacle
        SX distance = (pow(x_pos - obs(0), 2) + pow(y_pos - obs(1), 2));

        // Define the constraint (distance must be greater than d_safe)
        g_obs(k) = distance - d_safe;
    }

    // Initialize lower and upper bounds for the constraints
    o_lb = SX::zeros(N, 1); // g_obs should be non-negative
    o_ub = SX::inf(N, 1); // No upper limit

    for (size_t i = 0; i < params.N; i++)
    {
        args.lbg(i + n_states*(N+1)+N*4*4) = o_lb(i); //n_states*(N+1)+N*4*4+1
        args.ubg(i + n_states*(N+1)+N*4*4) = o_ub(i);
    }
    
}

SX h(const SX& st, const SX& obs, const SX& d_safe) {
    return mtimes(st(0) - obs(0), st(0) - obs(0)) + mtimes(st(1) - obs(1), st(1) - obs(1)) - mtimes(d_safe, d_safe);
}

// calc updates psi2, passed by reference
void calc_cbf(const SX& st0, const SX& st1, const SX& st2, const SX& P, SX& psi2) {
        SX alpha1 = 0.1, alpha2 = 0.05;
        SX d_safe = 0.5;

        int end_idx_P = P.numel();
        SX obs = vertcat(P(end_idx_P - 2 - 16 - 12), P(end_idx_P - 1 - 16 - 12)); // Obstacle position 
                                                                        // Added 16 to account for the dynamic CBF

        SX psi1 = h(st1, obs, d_safe) - h(st0, obs, d_safe) + mtimes(alpha1, h(st0, obs, d_safe));
        psi2 = h(st2, obs, d_safe) - h(st1, obs, d_safe) + mtimes(alpha1, h(st1, obs, d_safe)) + mtimes(alpha2, psi1);
}

void cbf_constraints(const casadi::SX& X, const casadi::SX& P, const Params& params, casadi::SX& g_obs, Args& args, const SX& slack_var) {
    using namespace casadi;

    casadi::SX o_lb, o_ub;
    int n_states = params.n_states;
    int N = params.N; 
    SX psi2; // this is the output of calc_cbf function

    // Initialize constraint array
    g_obs = SX::zeros(params.N, 1); // Use vertcat for a column vector

    // Loop over each time step
    for (int k = 1; k < N-2; ++k) {
        SX st0 = X(Slice(), k);       // X(:, k)
        SX st1 = X(Slice(), k+1);     // X(:, k+1)
        SX st2 = X(Slice(), k+2);     // X(:, k+2)

        calc_cbf(st0, st1, st2, P, psi2); 
        g_obs(k-1) = psi2; // + slack_var;
    }

    // manually updating the last 2 elements of g_obs, assume N = 15
    // CBF only computes N-2 elements, and collision constraints are N elements for now
    
    SX vel_ref = SX::zeros(1, 1);
    SX vel = SX::zeros(1, 1);
    // calculate average velocity
    for (size_t i = 1; i <= N; i++)
    {
        int index3 = 3 + i * 12;  // Corresponds to X(3, i)
        int index4 = 4 + i * 12;  // Corresponds to X(4, i)

        vel = vel +  casadi::SX::pow(X(3, i),2) + casadi::SX::pow(X(4, i),2) ;
        vel_ref = vel_ref + SX::pow(P(index3), 2) + SX::pow(P(index4), 2);
    }

    vel = vel / N;
    vel_ref = vel_ref / N;

    SX pos_err = SX::zeros(1, 1);
    // calculate position error
    for (size_t i = 1; i <= N; i++)
    {
        int index3 = (i) * (12);  // Corresponds to X(3, i)

        pos_err = pos_err + X(3, i)*P(index3);
    }

    vel = vel / N;
    vel_ref = vel_ref / N;

    writeMatrixToFile(vel, "tmp/vel.txt");
    writeMatrixToFile(vel_ref, "tmp/vel_ref.txt");
    
    // Direction matching constraint using squared norms
    

    g_obs(N-2) = pos_err;
    g_obs(N-1) = (vel_ref - vel); // + slack_var;

    // Initialize lower and upper bounds for the constraints
    o_lb = SX::zeros(N, 1); // g_obs should be non-negative
    o_ub = SX::inf(N, 1); // No upper limit

    for (size_t i = 0; i < params.N; i++)
    {
        args.lbg(i + n_states*(N+1)+N*4*4) = o_lb(i); //n_states*(N+1)+N*4*4+1
        args.ubg(i + n_states*(N+1)+N*4*4) = o_ub(i);
    }

    writeMatrixToFile(g_obs, "tmp/g_obs.txt");
    
}





#endif // NMPC_H