// Author: Basit Muhammad Imran 
// Date: 2024-04-17

#include "astar.hpp"
#include "nmpc.hpp"

NMPC::~NMPC() {
    xoutFile_.close();
    trajFile_.close();
    uoutFile_.close();
    feetFile_.close();
    toutFile_.close();
    cbfFile_.close();
    nmpcSolveTimeFile_.close();
}

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

casadi::DM NMPC::bezier(const DM& alpha, double s) {
    int n = alpha.size1(); // Number of rows
    int m = alpha.size2(); // Number of columns
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
            value(i) += alpha(i, j) * k[j] * x[j] * y[M - j];
        }
    }

    return value;
}

DM NMPC::bezierd(const DM& alpha, double s) {
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

SX NMPC::Rotz(const SX& t) {
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

casadi::SX NMPC::skewSym(const casadi::SX& x) {
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

void NMPC::planner() {
    const int n_states = params.n_states; // Example: Number of states
    // const int N = params.N;
    
    DM xd = DM::zeros(n_states * params.N);
    DM xd_log = DM::zeros(n_states * params.N);

    if (params.current_time <= 1.0) { // Stand up duration
        params.gait = 1;
    }
    else {
        params.gait = 8;
    }

    // writeMatrixToFile(q_des_, basePath + "q_des_.txt");

    DM traj = q_des_(casadi::Slice(0, n_states), 0); // Second parameter 0 is to indicate column index for a vector


    switch (params.gait) {
        case 1: {
            params.pf(0,0) = toePos_(0, 0) + agent_Initial_(0);  params.pf(1,0) = toePos_(1, 0) + agent_Initial_(0); params.pf(2,0) = toePos_(2, 0);
            params.pf(0,1) = toePos_(0, 1) + agent_Initial_(0);  params.pf(1,1) = toePos_(1, 1) + agent_Initial_(0); params.pf(2,1) = toePos_(2, 1);
            params.pf(0,2) = toePos_(0, 2) + agent_Initial_(0);  params.pf(1,2) = toePos_(1, 2) + agent_Initial_(0); params.pf(2,2) = toePos_(2, 2);
            params.pf(0,3) = toePos_(0, 3) + agent_Initial_(0);  params.pf(1,3) = toePos_(1, 3) + agent_Initial_(0); params.pf(2,3) = toePos_(2, 3);
            double t0 = 0, tf = 1.0; // Start and end times for mode 1
            double s = std::min((params.current_time - t0) / (tf - t0), 1.0); // Ensure s is within [0, 1]
            
            // Extract relevant parameters for easy access
            double xf = agent_Initial_(0); // Initial x position
            double yf = agent_Initial_(1); // Initial y position
            double x0 = xf;
            double y0 = yf;
            double zf = params.height; // Assuming the height is constant or predefined
            double z0 = 0.08;

            // Using DM::horzcat to create each row of the alpha matrix
            casadi::DM row1 = casadi::DM::horzcat({DM(x0), DM(x0), DM(x0), DM(xf / 4), DM(3 * xf / 4), DM(xf), DM(xf), DM(xf)});
            casadi::DM row2 = casadi::DM::horzcat({DM(y0), DM(y0), DM(y0), DM(yf / 4), DM(3 * yf / 4), DM(yf), DM(yf), DM(yf)});
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
            DM Rz = Rotz(vertcat(qk_(6), qk_(7), qk_(8))); // Assuming q contains Euler angles in the last 3 elements
            // size_t currentStep = std::floor(t / params.dt) + 1;
            // std::cout << "What is NMPC domain?: " << domain_ << std::endl;
            // std::cout << "What is params.stpcnt?: " << params.stepcnt << std::endl;
            // std::cout << "What is params.current_time: " << params.current_time << std::endl;
            // std::cout << "What is control_tick_: " << control_tick_*0.001 << std::endl;
            // std::cout << "What is my toePos_: " << toePos_ << std::endl;

            // TODO-Integration : LL is responsible for handling the domain change and corresponding updates
            // Such as changing the contact sequence, increasing velocity, and updating step foot position
            // if (params.current_time >= (params.stepcnt * domLenSec + 1)) {  // Domain change condition // params.stepcnt previously defined

            // TODO-Integration : Enforce domain change through LL
            //if (params.current_time >= (domain_ * domLenSec + 1)) {
            if (isDomainChanged_) {
                // Update contact sequence
                if (static_cast<double>(params.contact(0)) == 1.0) {
                    params.contact = {0, 1, 1, 0};
                } else {
                    params.contact = {1, 0, 0, 1};
                }

                // // Increase velocity in body frame
                // traj(3) += 0.1; // Increase command velocity in body X
                // traj(4) += 0.1; // Increase command velocity in body Y

                // // Saturate the command velocities
                // traj(3) = fmin(traj(3), params.desVel(0));
                // traj(4) = fmin(traj(4), params.desVel(1));

                // traj(3) = Prd_refined_(0, params.mpc_iter);
                // traj(4) = Prd_refined_(1, params.mpc_iter);


                // Calculate step length in body frame and update step foot position
                // DM stepLen =  vertcat(traj(3) * domLenSec / 2, 0.0, 0.0);;// vertcat(traj(3) * domLenSec / 2, traj(4) * domLenSec / 2, 0.0);
                
                
                // Calculate new-footholds
                params.pf(0,0) = toePos_(0, 0);  params.pf(1,0) = toePos_(1, 0); params.pf(2,0) = toePos_(2, 0);
                params.pf(0,1) = toePos_(0, 1);  params.pf(1,1) = toePos_(1, 1); params.pf(2,1) = toePos_(2, 1);
                params.pf(0,2) = toePos_(0, 2);  params.pf(1,2) = toePos_(1, 2); params.pf(2,2) = toePos_(2, 2);
                params.pf(0,3) = toePos_(0, 3);  params.pf(1,3) = toePos_(1, 3); params.pf(2,3) = toePos_(2, 3);

                // params.p_hip(0,0) = toePos_(0, 0);  params.p_hip(1,0) = toePos_(1, 0); params.p_hip(2,0) = toePos_(2, 0);
                // params.p_hip(0,1) = toePos_(0, 1);  params.p_hip(1,1) = toePos_(1, 1); params.p_hip(2,1) = toePos_(2, 1);
                // params.p_hip(0,2) = toePos_(0, 2);  params.p_hip(1,2) = toePos_(1, 2); params.p_hip(2,2) = toePos_(2, 2);
                // params.p_hip(0,3) = toePos_(0, 3);  params.p_hip(1,3) = toePos_(1, 3); params.p_hip(2,3) = toePos_(2, 3);

                // DM hipTmp = params.p_hip + stepLen;
                // hipTmp = mtimes(Rz, hipTmp); // Rotate to the world frame
                // params.pf = hipTmp + vertcat(qk_(Slice(0, 2)), 0.0);

                // Set velocity in world frame
                // DM cmdWrld = mtimes(Rz, traj(Slice(3, 6))); // Convert to world frame
                // traj(Slice(3, 6)) = cmdWrld;
                // params.stepcnt++;
                params.stepcnt = domain_;
            }

            int idx = findNearestIndex();

            // if (params.mpc_iter >= 100) {
            //     idx = params.mpc_iter;    
            // }
            // else {
            //     idx = 0;
            // }

            traj(0) =  Pr_refined_(agent_id_ * 2, idx);
            traj(1) =  Pr_refined_(agent_id_ * 2 + 1, idx);
            traj(3) = Prd_refined_(agent_id_ * 2, idx);
            traj(4) = Prd_refined_(agent_id_ * 2 + 1, idx);

            float heading = 0;

            // Update position and height based on velocity
            // traj(0) = qk_(0) + params.dt * traj(3);
            // traj(1) = qk_(1) + params.dt * traj(4);
            traj(2) = params.height; // Assuming height is constant or predefined
            traj(5) = 0; // Assuming zero vertical velocity

            for (int n = 1; n <= (int)params.N; ++n) {
                // Compute the predicted state for each step in the horizon
                DM state = DM::zeros(n_states, 1);
                DM state_log = DM::zeros(n_states, 1);

                // Index to access the right column in Prd_refined_
                int index = params.mpc_iter + n - 1;  // Adjust index for zero-based indexing if necessary

                // Make sure we don't exceed the bounds of Prd_refined_
                if (index >= Prd_refined_.cols()) {
                    std::cerr << "Index exceeds the size of Prd_refined_" << std::endl;
                    break;  // or handle overflow more gracefully
                }

                // Position components from Prd_refined_
                double cx = Pr_refined_(agent_id_ * 2, idx + 1 + n);  // Access the velocity in x-direction
                double cy = Pr_refined_(agent_id_ * 2 + 1, idx + 1 + n);  // Access the velocity in y-direction

                // TODO : add velocity component to converge back to the desired position
                Eigen::Vector2d pos_err = Eigen::Vector2d(qk_(0)-cx, qk_(1)-cy);

                double vx_corr = 0;
                double vy_corr = 0;

                // find min_obs
                Eigen::Vector2d obs_dist; obs_dist.setZero();
                obs_dist << (double)qk_(0) - obs_near_(0), (double)qk_(1) + obs_near_(1);

                // TODO : find min obstacle distance

                // FIXME : disable vel correction for now
                if (obs_dist.norm() > 0.75)
                {
                    vx_corr = ((double)qk_(0)-cx)*0.7;
                    vy_corr = ((double)qk_(1)-cy)*0.7;
                }

                std::cout << "pos_err: " << pos_err.transpose() << "    obs_dist:  " << obs_dist.norm() << std::endl;

                // Velocity components from Prd_refined_
                double vx = Prd_refined_(agent_id_ * 2, idx + 1 + n) - vx_corr;  // Access the velocity in x-direction
                double vy = Prd_refined_(agent_id_ * 2 + 1, idx + 1 + n) - vy_corr;  // Access the velocity in y-direction;

                // Updating the state based on velocities from Prd_refined_
                state(0) = qk_(0) + n * params.dt * vx;  // pos X
                state(1) = qk_(1) + n * params.dt * vy;  // pos Y

                state(2) = params.height;                // pos Z (constant height assumption)
                state(3) = vx;                           // vel X (updated from Prd_refined_)
                state(4) = vy;                           // vel Y (updated from Prd_refined_)
                state(5) = 0.0;                          // vel Z (assuming zero vertical velocity)
                // Assuming Euler angles (heading) and angular velocity are zero
                state(6) = 0.0;  // Roll or phi
                state(7) = 0.0;  // Pitch or theta
                state(8) = heading;  // Yaw or psi
                state(9) = 0.0;      // Angular velocities (assuming zeros)
                state(10) = 0.0;
                state(11) = 0.0;

                state_log = state;
                // state_log(0) = cx;
                // state_log(1) = cy;

                // Place the computed state into the appropriate segment of xd
                // Adjust indexing if your container xd is one-dimensional but structured like a two-dimensional array
                xd(Slice((n-1)*n_states, n*n_states)) = state;
                // xd_log(Slice((n-1)*n_states, n*n_states)) = state_log;
            }

            break;
        }
        // Handle other cases as necessary
        default:
            std::cerr << "Mode '" << params.gait << "' is not defined in the motion planner.\n";
    }

    if (params.gait != 8) {
        xd = repmat(traj, params.N, 1);
        // xd_log = repmat(traj, params.N, 1);
    }

    q_des_ = xd; // Only xd needs to be returned
    // q_des_log_ = xd_log;
}
void NMPC::force_inequality(const casadi::SX& U, casadi::SX& g_force) {
    // Use SX for symbolic constraint expressions
    casadi::SX g_force_sx = casadi::SX::vertcat(std::vector<casadi::SX>{});

    // g_ub and g_lb as symbolic vectors with appropriate bounds
    casadi::DM g_ub_DM = casadi::DM::vertcat(std::vector<casadi::DM>{});
    casadi::DM g_lb_DM = casadi::DM::vertcat(std::vector<casadi::DM>{});


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
            g_ub_DM = casadi::SX::vertcat({g_ub_DM, 0.0, 0.0, casadi::inf, casadi::inf}); 
            g_lb_DM = casadi::SX::vertcat({g_lb_DM, -casadi::inf, -casadi::inf, 0, 0});
        }
    }


    // Convert symbolic SX expressions to DM for g_ub and g_lb if necessary
    // Note: This might not be needed if keeping them symbolic is preferable
    g_force = g_force_sx; // g_force is already SX and can stay symbolic
    g_ub_ = g_ub_DM;
    g_lb_ = g_lb_DM;
}

void NMPC::force_inequality_bounds() {

    // g_ub and g_lb as DM vectors with appropriate bounds
    casadi::DM g_ub_sx;
    casadi::DM g_lb_sx;

    for (size_t k = 0; k < params.N; ++k) {
        for (size_t leg = 0; leg < 4; ++leg) {
            // Directly construct and append symbolic expressions
            g_ub_sx = casadi::SX::vertcat({g_ub_sx, 0.0, 0.0, casadi::inf, casadi::inf}); 
            g_lb_sx = casadi::SX::vertcat({g_lb_sx, -casadi::inf, -casadi::inf, 0, 0});
        }
    }
    g_ub_ = g_ub_sx;
    g_lb_ = g_lb_sx;
}

void NMPC::writeMatrixToFile(const casadi::SX& matrix, const std::string& filename) {
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

void NMPC::appendMatrixToFile(const casadi::DM& matrix, std::ofstream& file) {
    if (file.is_open()) {
        for (int i = 0; i < matrix.size1(); ++i) {
            for (int j = 0; j < matrix.size2(); ++j) {
                file << matrix(i, j);
                if (j < matrix.size2() - 1) file << ", ";
            }
            file << std::endl;
        }
    } else {
        std::cerr << "Unable to write to file" << std::endl;
    }
}

void NMPC::writeMatrixToFile(const casadi::DM& matrix, const std::string& filename) {
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

void NMPC::writeMatrixToFile(const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::trunc);  // Open for output and truncate to clear existing content
    if (file.is_open()) {
        // Helper lambda to write a single matrix
        auto writeMatrix = [&](const casadi::DM& matrix) {
            // file << "Matrix size: " << matrix.size1() << "x" << matrix.size2() << "\n" << std::endl;
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

void NMPC::equalityconstraints(const SX& X, const SX& U, const SX& P, SX& obj, SX& g_eq, const Function& f) {

    std::string basePath = "tmp/";
    int n_states = params.n_states;
    double dt = params.dt;

    // Initial condition constraints
    SX st = X(Slice(0, n_states), 0);
    SX pf = P(Slice(P.numel()-12, 12));
    SX footholds = reshape(pf, 4, 3); 
    g_eq = st - P(Slice(0, n_states));

    for (int k = 0; k < (int)params.N; ++k) {
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
        SX stage_cost = mtimes(mtimes(st_err.T(), params.Q), st_err);
        SX terminal_cost = mtimes(mtimes(con.T(), params.R), con);

        // writeMatrixToFile(stage_cost, basePath + "stage_cost.txt");
        // writeMatrixToFile(terminal_cost, basePath + "terminal_cost.txt");

        obj = obj + stage_cost + terminal_cost;

        // Calculate next state prediction
        SX st_next = X(Slice(0, n_states), k + 1);

        std::map<std::string, SX> arg;
        arg["st"] = st;  // Example random inputs
        arg["con"] = con;
        arg["pf"] = footholds;

        std::map<std::string, SX> result = f(arg);
        SX f_value = result.at("rhs");

        // writeMatrixToFile(st_next, basePath + "st_next.txt");
        // writeMatrixToFile(f_value, basePath + "f_value.txt");

        // Compute constraints
        SX st_next_euler = st + dt * f_value;
        g_eq = SX::vertcat({g_eq, st_next - st_next_euler});

        // writeMatrixToFile(st_next_euler, basePath + "st_next_euler.txt");
        // writeMatrixToFile(g_eq, basePath + "g_eq.txt");
    }
}

void NMPC::set_constraint_bounds() {
    // Calculate the start index for setting lower and upper bounds

    int startIndex = params.n_states * (params.N + 1); // Adjust for zero-based indexing

    // Resize the lbg and ubg if necessary
    if (args.lbg.is_empty() || args.lbg.size1() < startIndex + g_lb_.size1()) {
        args.lbg = DM::zeros(startIndex + g_lb_.size1(), 1); // Adjust total size if needed
    }
    if (args.ubg.is_empty() || args.ubg.size1() < startIndex + g_ub_.size1()) {
        args.ubg = DM::zeros(startIndex + g_ub_.size1(), 1); // Adjust total size if needed
    }

    // Set values for lower bounds
    for (int i = 0; i < g_lb_.size1(); ++i) {
        args.lbg(startIndex + i) = g_lb_(i);
    }

    // Set values for upper bounds
    for (int i = 0; i < g_ub_.size1(); ++i) {
        args.ubg(startIndex + i) = g_ub_(i);
    }
}

void NMPC::set_statebounds() {
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
        for (int leg = 0; leg < 4; ++leg) {
            int index = n_states*(N+1) + n_inputs*k + 3*leg + 2; // Adjusted for 0-based indexing
            args.lbx(index) = 0;
            args.ubx(index) = params.contact(leg) * zmax; // Assuming contact is std::vector<double>
        }
    }
}

void NMPC::shift() {
    // Prepare input for the function f
    std::map<std::string, DM> arg;
    arg["st"] = qk_;  // Current state
    arg["con"] = ukAll_(0, Slice()).T();  // First control input, transposed to match MATLAB

    // Call the function f
    std::map<std::string, DM> result = f_dyn_(arg);
    DM f_value = result["rhs"];  // Assuming the output is named 'rhs'

    // Update state
    // qk_ = qk_ + DM(params.dt) * f_value;

    // Increment time
    // params.current_time += params.dt;

    DM prev_ukAll = ukAll_;

    // Shift control inputs
    DM u_new = vertcat(prev_ukAll(Slice(1, prev_ukAll.size1()), Slice()),  // Skip the first control set
                       prev_ukAll(Slice(prev_ukAll.size1()-1, prev_ukAll.size1()), Slice()));  // Duplicate the last control set
    ukAll_ = u_new;
}

void NMPC::code_gen() {

    const double dt = params.dt; // Sample time
    const double sim_time = params.sim_length; // Simulation time
    const float sim_length = sim_time/dt; // Simulation length
    const int n_states = params.n_states; // Number of states
    const int n_inputs = params.n_inputs; // Number of inputs
    const int n_controls = params.n_controls; // Number of inputs
    const int gait = params.gait; // Gait type
    const double G = params.G; // Gravity
    const double mass = params.mass; // Mass

    init(); // Initialize the variables

    // Symbolic variables 
    casadi::SXVector states_vec = {
        SX::sym("Xx"), SX::sym("y"), SX::sym("z"),
        SX::sym("dx"), SX::sym("dy"), SX::sym("dz"),
        SX::sym("phi"), SX::sym("theta"), SX::sym("psi"),
        SX::sym("dphi"), SX::sym("dtheta"), SX::sym("dpsi")
    };

    // Define controls similarly
    casadi::SXVector controls_vec = {
        SX::sym("u1x"), SX::sym("u1y"), SX::sym("u1z"),
        SX::sym("u2x"), SX::sym("u2y"), SX::sym("u2z"),
        SX::sym("u3x"), SX::sym("u3y"), SX::sym("u3z"),
        SX::sym("u4x"), SX::sym("u4y"), SX::sym("u4z")
    };

    casadi::SX states = casadi::SX::vertcat(states_vec);
    casadi::SX controls = casadi::SX::vertcat(controls_vec);

    SX U = SX::sym("U", n_controls, params.N);
    SX P = SX::sym("P", n_states + n_states * params.N);
    SX X = SX::sym("X", n_states, params.N + 1);
    SX OPT_variables = vertcat(reshape(X, n_states * (params.N + 1), 1), reshape(U, n_inputs * params.N, 1));


    DM qk = DM::zeros(n_states);
    qk(2) = 0.08; // Assuming z position is constant or predefined

    double t0 = 0.0;
    DM xx(n_states, sim_length);
    xx(Slice(0, n_states), 0) = qk;  

    DM u0 = DM::zeros(params.N, params.n_inputs); // Control inputs for N horizon steps, not sim_steps

    // States decision variables initialization
    DM X0 = repmat(qk, 1, params.N+1);
    DM xd = repmat(qk, params.N, 1);

    double current_time = 0.0;
    size_t mpciter = 0;

    Args args;
    args.lbg = DM::zeros(n_states * (params.N + 1) + 4*4*params.N, 1);
    args.ubg = DM::zeros(n_states * (params.N + 1) + 4*4*params.N, 1);

    DM uout = DM::zeros(n_controls, sim_length); // Vector to store control outputs
    DM xout = DM::zeros(n_states, sim_length);
    DM feet = DM::zeros(12, sim_length); // Vector to store foot positions
    DM traj = DM::zeros(12, sim_length); // Vector to store desired trajectory
    DM tout = DM::zeros(1, sim_length);
    DM cbfout = DM::zeros(1, sim_length);


    // params.current_time = current_time;
    params.current_time = control_tick_*0.001;
    params.mpc_iter = mpciter;
    planner(); // Assuming gait and other details are handled inside
    args.p = vertcat(qk, xd);

    SX t_rot = casadi::SX::vertcat({states(6), states(7), states(8)});  // Extract phi, theta, and psi into a column vector for Rotz
    SX Rz = Rotz(t_rot);    // Call Rotz with the column vector
    SX z3 = casadi::SX::zeros(3, 3);    // Create a 3x3 zero matrix
    SX i3 = casadi::SX::eye(3);         // Create a 3x3 identity matrix

    DM J = DM::zeros(3,3); // Assuming a paramse-5; J(2,2) = 0.064713601;
    SX omega = casadi::SX::vertcat({states(9), states(10), states(11)}); // Extracting angular velocities (dphi, dtheta, dpsi)
    SX skew_omega = skewSym(omega);
    SX J_sx = SX(J); // Convert J from DM to SX
    SX IwI = mtimes(mtimes(mtimes(Rz, J_sx), Rz.T()), skew_omega) * mtimes(mtimes(Rz, J_sx), Rz.T());

    SX A = SX::zeros(12, 12); // A is a 12x12 matrix given your block structure
    // Fill in the blocks of A
    A(Slice(0, 3), Slice(0, 3)) = z3;
    A(Slice(0, 3), Slice(3, 6)) = i3;
    A(Slice(0, 3), Slice(6, 9)) = z3;
    A(Slice(0, 3), Slice(9, 12)) = z3;

    A(Slice(3, 6), Slice(0, 3)) = z3;
    A(Slice(3, 6), Slice(3, 6)) = z3;
    A(Slice(3, 6), Slice(6, 9)) = z3;
    A(Slice(3, 6), Slice(9, 12)) = z3;

    A(Slice(6, 9), Slice(0, 3)) = z3;
    A(Slice(6, 9), Slice(3, 6)) = z3;
    A(Slice(6, 9), Slice(6, 9)) = z3;
    A(Slice(6, 9), Slice(9, 12)) = Rz.T();

    A(Slice(9, 12), Slice(0, 3)) = z3;
    A(Slice(9, 12), Slice(3, 6)) = z3;
    A(Slice(9, 12), Slice(6, 9)) = z3;
    A(Slice(9, 12), Slice(9, 12)) = z3; //-IwI
    
    SX D = SX::zeros(12,1); // D is a 12x1 vector
    D(Slice(3,6)) = SX::vertcat({SX::zeros(2,1), -G}); // Fill the 4th to 6th rows

    SX p_foot_sx = SX(params.pf);
    SX com_vector = SX::vertcat({states(0), states(1), states(2)});
    SX com_vector_expanded = repmat(com_vector, 1, 4); // Replicate across columns to match p_foot_sx dimensions

    SX rd = p_foot_sx - com_vector_expanded;
    
    SX skew_rd = skewSym(rd);
    SX B_tmp = mtimes(mtimes(Rz, J_sx), Rz.T()) * skew_rd;

    SX B = SX::vertcat({
        SX::zeros(3, 12),
        SX::repmat(1/mass * i3, 1, 4),
        SX::zeros(3, 12),
        mtimes(mtimes(Rz, J_sx), Rz.T()) * skew_rd}); //mtimes(mtimes(Rz, J_sx), Rz.T()) * skew_rd
    SX rhs = mtimes(A, states) + mtimes(B, controls) + D;

    // writeMatrixToFile(rhs, basePath + "rhs.txt");

    Function f_dyn = Function("f_dyn", {states, controls}, {rhs}, {"st", "con"}, {"rhs"});

    // Declare variables to hold the outputs
    SX g_force;

    // Call the function
    force_inequality(U, g_force);

    SX obj = 0, g_eq; 
    equalityconstraints(X, U, P, obj, g_eq, f_dyn);
    SX g = vertcat(g_eq, g_force);
    set_statebounds();
    set_constraint_bounds();
    SXDict nlp_prob = {{"f", obj}, {"x", OPT_variables}, {"g", g}, {"p", P}};

    // Define options for the solver
    Dict opts;

    // For ipopt
    opts["ipopt.max_iter"] = 10;  // Replace Max_mpciter with its actual value
    opts["ipopt.print_level"] = 0;  // Can be changed to 0 or 3 based on verbosity required
    opts["print_time"] = 0;  // Disable printing solver time
    opts["ipopt.acceptable_tol"] = 1e-2;  // Tolerance for stopping criterion
    opts["ipopt.acceptable_obj_change_tol"] = 1e-2;  // Objective change tolerance for stopping

    /// for snopt
    // opts["snopt.Iterations limit"] = 10;
    // opts["snopt.Major Print level"] = 3;
    // opts["snopt.Minor Print level"] = 3;
    // opts["snopt.Derivative option"] = 0;
    // opts["snopt.Major feasibility tolerance"] = 1.0e-3;
    // opts["snopt.Minor feasibility tolerance"] = 1.0e-3;
    // opts["snopt.Major optimality tolerance"] = 1.0e-3;
    // opts["snopt."] = 10;
    // opts["print_level"] = 0;
    // // opts["print_file"] = "snopt.out";
    // opts["major_feasibility_tolerance"] = 1e-2;
    // opts["major_optimality_tolerance"] = 1e-2;

    // Create an NLP solver instance
    Function solver = nlpsol("solver", "ipopt", nlp_prob, opts);
    // Function solver = nlpsol("solver", "snopt", {{"f", obj}, {"x", OPT_variables}, {"g", g}, {"p", P}});

    // Reshape X0 and u0 transposed into column vectors
    DM X0_col = reshape(X0.T(), n_states * (params.N + 1), 1);
    DM u0_col = reshape(u0.T(), n_inputs * params.N, 1);

    // Concatenate reshaped vectors vertically
    DM x0u0= vertcat(X0_col, u0_col);
    args.x0 = x0u0;
    args.p = vertcat(qk, xd);
    // Solve the problem
    DMDict arg = {{"x0", args.x0}, {"lbx", args.lbx}, {"ubx", args.ubx}, {"lbg", args.lbg}, {"ubg", args.ubg}, {
        "p", args.p}};
    DMDict res = solver(arg);


    DM sol = DM::zeros(n_states * (params.N + 1) + n_controls * params.N, 1); // sol.x should be initialized with actual solution data
    sol = res.at("x");


    // =========================== CODE GEN =========================== //

    std::string file_name = "nlp_code";
    // code predix
    std::string prefix_code = fs::current_path().string() + "/";
    casadi::Dict code_gen_opts = casadi::Dict();
    code_gen_opts["cpp"] = false;
    code_gen_opts["with_header"] = false;
    solver.generate_dependencies(file_name + ".c");

    CodeGenerator codegen = CodeGenerator(file_name + "_f" + ".c", code_gen_opts);
    codegen.add(f_dyn);
    codegen.generate();
    std::string prefix_lib = fs::current_path().string() + "/";
    std::string compile_command_solver = "gcc -fPIC -shared -O3 " + 
        prefix_code + file_name + ".c -o " +
        prefix_lib + file_name + ".so -L/usr/local/lib -lipopt";
    
    std::string compile_command_f = "gcc -fPIC -shared -O3 " + 
        prefix_code + file_name + "_f" + ".c -o " +
        prefix_lib + file_name + "_f" + ".so";
    std::cout << compile_command_solver << std::endl;
    std::cout << compile_command_f << std::endl;

    int compile_flag_solver = std::system(compile_command_solver.c_str());
    casadi_assert(compile_flag_solver==0, "Compilation failed");
    std::cout << "Compilation successed!" << std::endl;

    int compile_flag_f = std::system(compile_command_f.c_str());
    casadi_assert(compile_flag_f==0, "Compilation failed");
    std::cout << "Compilation successed!" << std::endl;
}

void NMPC::init()
{
    qk_ = DM::zeros(params.n_states);
    qk_(2) = 0.08; // sitting com height
    q_des_ = DM::zeros(params.n_states*params.N, 1);
    ukAll_ = DM::zeros(params.N, params.n_inputs);

    qkAll_ = repmat(qk_, 1, params.N+1); //contains initial state, and solution for all N horizon steps
    q_des_ = repmat(qk_, params.N, 1);   //contains desired state for all N horizon steps

    args.lbg = DM::zeros(params.n_states * (params.N + 1) + 4*4*params.N + params.N, 1);
    args.ubg = DM::zeros(params.n_states * (params.N + 1) + 4*4*params.N + params.N, 1);

    uout = DM::zeros(params.n_controls, 1);//params.sim_length*1000
    xout =  DM::zeros(params.n_states, 1);
    feet =  DM::zeros(12, 1);   // Vector to store foot positions   
    traj =  DM::zeros(12, 1); 
    tout =  DM::zeros(1, 1);
    cbfout =  DM::zeros(1, 1);
    NMPC_solve_time = DM::zeros(1, 1);  

    std::string file_name = "nlp_code";
    std::string prefix_lib = fs::current_path().string() + "/";
    
    std::string lib_name = prefix_lib + file_name + ".so";
    // Load the nlp solver
    // Define options for the solver-opts
    Dict opts;

    int iter_lim = 10;
    if (params.current_time < 1)
    {
        iter_lim = 8; // 13
    }
    else 
        iter_lim = 8; // 8, 12 also works
    
    opts["ipopt.max_iter"] = iter_lim;  // Replace Max_mpciter with its actual value
    opts["ipopt.print_level"] = 0;  // Can be changed to 0 or 3 based on verbosity required
    opts["print_time"] = 0;  // Disable printing solver time
    opts["ipopt.acceptable_tol"] = 1e-3;  // Tolerance for stopping criterion
    opts["ipopt.acceptable_obj_change_tol"] = 1e-3;  // Objective change tolerance for stopping
    solver_ = casadi::nlpsol("solver","ipopt", lib_name, opts);

    // Load nonlinear dynamics function 
    lib_name = prefix_lib + file_name + "_f" ".so";
    f_dyn_ = external("f_dyn", lib_name);
}

void NMPC::set_collision_bounds() {
    using namespace casadi;

    casadi::DM o_lb, o_ub;
    int n_states = params.n_states;
    int N = params.N;

    // Initialize lower and upper bounds for the constraints
    o_lb = DM::zeros(N, 1); // g_obs should be non-negative
    o_ub = DM::inf(N, 1); // No upper limit

    // for vel constraint
    // o_ub(N-2) = 0.15;
    o_ub(N-1) = 0.025;

    for (size_t i = 0; i < params.N; i++)
    {
        args.lbg(i + n_states*(N+1)+N*4*4) = o_lb(i); //n_states*(N+1)+N*4*4+1
        args.ubg(i + n_states*(N+1)+N*4*4) = o_ub(i);
    }
    
}

void NMPC::run_NMPC() {

    // call the planner function to update q_des_ for N horizon steps
    planner();

    // writeMatrixToFile(q_des_, basePath + "q_des_.txt");
    // writeMatrixToFile(basePath + "params.txt");

    force_inequality_bounds();     // Updates g_ub_ and g_lb_ with appropriate bounds
    set_statebounds();                     // Updates args.lbx and args.ubx with appropriate bounds
    set_constraint_bounds(); // Updates args.lbg and args.ubg with g_lb_ and g_ub_
    set_collision_bounds(); // Updates args.lbg and args.ubg with o_lb and o_ub

    // collision_constraints(const casadi::SX& X, casadi::SX& g_obs, Args& args);
    // writeMatrixToFile(g_ub, basePath + "g_ub.txt");
    // writeMatrixToFile(g_lb, basePath + "g_lb.txt");
    // writeMatrixToFile(params, basePath + "params.txt");

    solve_nlp(); // Solve the NMPC problem
    logData();
}

void NMPC::solve_nlp() {

    DM q0_col = reshape(qkAll_.T(), params.n_states * (params.N + 1), 1);
    DM u0_col = reshape(ukAll_.T(), params.n_inputs * params.N, 1);

    // Concatenate reshaped vectors vertically
    DM x0u0= vertcat(q0_col, u0_col, DM(0)); // slack_var = 0
    args.x0 = x0u0;
    find_nearest_obs(); 

    // find min_obs
    // obs_near_(0) = state_other_(0); obs_near_(1) = state_other_(1); // Hardcoded for now
    std::cout << "Nearest obstacle : (" << agent_id_ << "): " << obs_near_.transpose() << std::endl;

    // args.p = vertcat(qk_, q_des_, DM(5), DM(-0.1)); // Add obs info here 
    DM other_agent_info = DM::zeros(16, 1);

    // DEFINE or PROPERLY CHANNEL FOOTHOLDS
    DM footholds = DM::zeros(12, 1);

    // Convert Eigen::VectorXd to casadi::DM
    casadi::DM GRF_sol_other_DM = casadi::DM::zeros(12, 1);
    casadi::DM state_other_DM = casadi::DM::zeros(4, 1);
    eigenToCasadi(GRF_sol_other_, GRF_sol_other_DM);
    eigenToCasadi(state_other_, state_other_DM);

    // std::cout << "state_other_: " << state_other_.transpose() << std::endl;
    // std::cout << "state_other_DM: " << state_other_DM << std::endl;

    // std::cout << "GRF_sol_other_: " << GRF_sol_other_DM << std::endl;
    // std::cout << "GRF_sol_other_DM: " << GRF_sol_other_DM << std::endl;


    if (GRF_sol_other_DM.numel() != 12) { // Check if the number of elements is 12
        // other_agent_info = casadi::DM::vertcat({GRF_sol_other_DM, state_other_DM});
        GRF_sol_other_DM = casadi::DM::zeros(12, 1);
    }

    other_agent_info = casadi::DM::vertcat({GRF_sol_other_DM, state_other_DM});

    // std::cout << "Other agent info: " << other_agent_info << std::endl;
    // std::cin.get();

    args.p = vertcat(qk_, q_des_, obs_near_(0), obs_near_(1), other_agent_info, footholds);

    // std::cout << "other_agent_info: " << args.p << std::endl;

    std::map<std::string, casadi::DM> arg, res;
    arg["lbx"] = args.lbx;
    arg["ubx"] =  args.ubx;
    arg["lbg"] =  args.lbg;
    arg["ubg"] =  args.ubg;
    arg["x0"] = args.x0;
    arg["p"] = args.p;

    // writeMatrixToFile(args.lbx, basePath + "args_lbx.txt");
    // writeMatrixToFile(args.ubx, basePath + "args_ubx.txt");
    // writeMatrixToFile(args.lbg, basePath + "args_lbg.txt");
    // writeMatrixToFile(args.ubg, basePath + "args_ubg.txt");
    // writeMatrixToFile(args.x0, basePath + "args_x0.txt");
    // writeMatrixToFile(args.p, basePath + "args_p.txt");
    
    auto start = std::chrono::high_resolution_clock::now();
    res = solver_(arg);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);


    DM sol = DM::zeros(params.n_states * (params.N + 1) + params.n_controls * params.N, 1); // sol.x should be initialized with actual solution data
    sol = res.at("x");

    DM g = res.at("g");
    
    // writeMatrixToFile(g, basePath + "g.txt");
    // writeMatrixToFile(args.p, basePath + "args_p.txt");

    casadi_int start_index_u = static_cast<casadi_int>(params.n_states * (params.N + 1));
    casadi_int end_index = static_cast<casadi_int>(sol.numel()-1);

    DM ukAll_sol = reshape(sol(Slice(start_index_u, end_index)).T(), params.n_controls, params.N).T();
    DM sol_x_N = reshape(sol(Slice(0, static_cast<casadi_int>(params.n_states * (params.N + 1)))).T(), params.n_states, params.N + 1).T();
    // Assuming sol_x_N is a previously defined CasADi DM or MX matrix
    DM mpc_sol_wo_IC = sol_x_N(Slice(1, sol_x_N.size1()), Slice());  // Slices all columns
    DM mpc_sol_w_IC = sol_x_N(Slice(0, sol_x_N.size1()), Slice());  // Slices all columns

    x_sol_N_ = mpc_sol_wo_IC; // assign class variable
    x_sol_N_w_IC = mpc_sol_w_IC; // assign class variable

    // writeMatrixToFile(x_sol_N_, basePath + "x_sol_N_.txt");
    // writeMatrixToFile(x_sol_N_w_IC, basePath + "x_sol_N_w_IC.txt");

    ukAll_ = ukAll_sol;


    // writeMatrixToFile(ukAll_, basePath + "ukAll_.txt");

    // writeMatrixToFile(ukAll_sol, basePath + "ukAll_.txt");
    // writeMatrixToFile(sol_x_N, basePath + "sol_x_N.txt");
    params.mpc_iter++;

    DM uk = ukAll_sol(0, Slice());
    uk_ = uk;
    uout(Slice(), 0) = uk;     // Store current control outputs
    
    // Enable the state-feedback here
    shift();
    qkAll_ = vertcat(sol_x_N(Slice(1, sol_x_N.size1()), Slice()), sol_x_N(Slice(sol_x_N.size1() - 1, sol_x_N.size1()), Slice()));
    // qk_ = qk_ + DM(params.dt) * f_value;
    
    // q:      [x, y, z, phi, theta, psi, ...]
    // dq:     [dx, dy, dz, dphi, dtheta, dpsi, ...]
    // qk_:    [x, y, z, dx, dy, dz, phi, theta, psi, dphi, dtheta, dpsi]

    // Positions and Velocities
    for (size_t i = 0; i < 3; i++) {
        qk_(i) = q[i];       // q: x, y, z
        qk_(i+3) = dq[i];    // dq: dx, dy, dz
    }

    qk_(2) = 0.28; qk_(5) = 0;
    // qk_(1) = q[1] - 0.005;

    // Rotation and Angular Velocities
    for (size_t i = 0; i < 3; i++) {
        qk_(i+6) = q[i+3];      // q: phi, theta, psi
        qk_(i+9) = dq[i+3];     // dq: dphi, dtheta, dpsi
    }
    
    

    DM sol_vector = reshape(x_sol_N_.T(), x_sol_N_.T().size1() * x_sol_N_.T().size2(), 1);

    // qkAll_ = reshape(sol_x_N, params.N+1*params.n_states, 1);
    x_sol_ = sol_vector;

    // writeMatrixToFile(x_sol_, basePath + "x_sol_.txt");

    // Data Logging
    // std::cout << "Starting Data Recording...\n";
    std::cout << "MPC Iter: " << params.mpc_iter << std::endl;
    xout(Slice(), 0) = qk_.T(); // Store current states, transposed
    traj(Slice(), 0) = q_des_(casadi::Slice(0, 12), 0).T(); // Store reshaped desired trajectory
    feet(Slice(), 0) = reshape(params.pf, 1, 12); // Store reshaped foot positions
    tout(0) = params.current_time - params.dt;          // Store current time
    NMPC_solve_time(0) = DM(duration.count()); // Store NMPC solve time
    cbfout(0) = (pow( qk_(0) - obs_near_(0), 2) + pow( qk_(1) - obs_near_(1), 2) - 0.65*0.65); // Store distance to nearest obstacle

    // std::cout << "Entering logData()...\n";
    logData();
    // std::cout << "Exiting logData()...\n";
    if (agent_id_ == 1) {
        // writeMatrixToFile(x_sol_N_, basePath + "x_sol_N_.txt");
        // writeMatrixToFile(traj, basePath + "traj.txt");
    }
    // writeMatrixToFile(ukAll_sol, basePath + "ukAll_sol.txt");
    // writeMatrixToFile(uk.T(), basePath + "uk_.txt");

    fitComTrajectory();
    fitGRFTrajectory();

    casadiToEigen(x_sol_, x_sol_eigen_);
    casadiToEigen(uk, GRF_sol_eigen_);
    casadiToEigen(bez_alpha_, bez_alpha_eigen_);
    // casadiToEigen(x_sol_N_w_IC.T(), bez_alpha_eigen_);
    casadiToEigen(bez_GRF_, bez_GRF_eigen_);

    // writeMatrixToFile(bez_GRF_, basePath + "bez_GRF_.txt");
    
    std::cout << "Current time: " << params.current_time << std::endl;
    params.current_time = control_tick_*0.001;


    // params.current_time += params.dt; is already done in shift()
}

// void NMPC::logData()
// {
//     std::string basePath = "tmp/";
//     writeMatrixToFile(xout.T(), basePath                   + "xout_" + std::to_string(agent_id_) + ".txt");
//     writeMatrixToFile(traj.T(), basePath                   + "traj_" + std::to_string(agent_id_) + ".txt");
//     writeMatrixToFile(uout.T(), basePath                   + "uout_" + std::to_string(agent_id_) + ".txt");
//     writeMatrixToFile(feet.T(), basePath                   + "feet_" + std::to_string(agent_id_) + ".txt");
//     writeMatrixToFile(tout.T(), basePath                   + "tout_" + std::to_string(agent_id_) + ".txt");
//     writeMatrixToFile(cbfout.T(), basePath                  + "cbf_" + std::to_string(agent_id_) + ".txt");
//     writeMatrixToFile(NMPC_solve_time, basePath + "NMPC_solve_time_" + std::to_string(agent_id_) + ".txt");
// }

void NMPC::logData() {
    appendMatrixToFile(xout.T(), xoutFile_);
    appendMatrixToFile(uout.T(), uoutFile_);
    appendMatrixToFile(traj.T(), trajFile_);
    appendMatrixToFile(feet.T(), feetFile_);
    appendMatrixToFile(tout, toutFile_);
    appendMatrixToFile(cbfout, cbfFile_);
    appendMatrixToFile(NMPC_solve_time, nmpcSolveTimeFile_);
}

int NMPC::binomial(int n, int k) {
        int result = 1;
        if (k > n - k) k = n - k;
        for (int i = 0; i < k; ++i) {
            result *= (n - i);
            result /= (i + 1);
        }
    return result;
}

int NMPC::binomialCoeff(int n, int k) {
        if (k > n) return 0;
        if (k == 0 || k == n) return 1;
        return binomialCoeff(n - 1, k - 1) + binomialCoeff(n - 1, k);
}

void NMPC::fitComTrajectory() {

    int degree = 4;  // Degree of the Bezier curve
    int num_points = params.N;
    int num_states = params.n_states;

    // Initialize matrix A for Bezier basis coefficients
    DM A = DM::zeros(num_points, degree + 1);

    // Function to calculate binomial coefficient
    auto binomialCoeff = [](int n, int k) -> double {
        double res = 1;
        if (k > n - k) k = n - k;
        for (int i = 0; i < k; ++i) {
            res *= (n - i);
            res /= (i + 1);
        }
        return res;
    };

    // Populate matrix A with Bezier basis values
    for (int i = 0; i < num_points; ++i) {
        double t = static_cast<double>(i) / (num_points - 1);  // Normalized time for this index
        for (int j = 0; j <= degree; ++j) {
            A(i, j) = binomialCoeff(degree, j) * pow(t, j) * pow(1 - t, degree - j);
        }
    }

    // Bezier coefficients matrix (one row per state)
    DM bez_alpha = DM::zeros(num_states, degree + 1);

    // Check if x_sol_N_ is empty or has incompatible dimensions
    if (x_sol_N_.is_empty() || x_sol_N_.rows() != num_points || x_sol_N_.columns() != num_states) {
        std::cerr << "Error: x_sol_N_ is empty or has incompatible dimensions." << std::endl;
        return;
    }

    // Modify matrix A and x_sol_N_ to account for the initial condition
    DM A_modified = A;
    DM x_sol_modified = x_sol_N_;
    
    // Modify the first row of A to be all zeros except for the first element
    A_modified(0, 0) = 1.0;
    for (int j = 1; j <= degree; ++j) {
        A_modified(0, j) = 0.0;
    }
    
    // Ensure the first row of x_sol_N_ starts from qk_
    for (int state = 0; state < num_states; ++state) {
        x_sol_modified(0, state) = qk_(state);
    }

    // Perform matrix operations for all states
    DM ATA = mtimes(transpose(A_modified), A_modified);
    DM ATx = mtimes(transpose(A_modified), x_sol_modified);  // Ensure correct dimension alignment

    bez_alpha = mtimes(pinv(ATA), ATx);

    bez_alpha_ = bez_alpha.T();

    // Output the Bezier Coefficients Matrix
    // writeMatrixToFile(bez_alpha_, basePath + "bez_alpha.txt");
    // std::cout << "Bezier Coefficients Matrix:" << std::endl;
    // std::cout << bez_alpha_ << std::endl;
    // std::cin.get();
}

void NMPC::fitGRFTrajectory() {
    int degree = 4;  // Degree of the Bezier curve
    int num_points = params.N;
    int num_states = params.n_states; // should be n_controls

    // Initialize matrix A for Bezier basis coefficients
    DM A = DM::zeros(num_points, degree + 1);

    // Populate matrix A with Bezier basis values
    for (int i = 0; i < num_points; ++i) {
        double t = static_cast<double>(i) / (num_points - 1);  // Normalized time for this index
        for (int j = 0; j <= degree; ++j) {
            A(i, j) = binomialCoeff(degree, j) * pow(t, j) * pow(1 - t, degree - j);
        }
    }

    // Bezier coefficients matrix (one row per state)
    DM bez_GRF = DM::zeros(num_states, degree + 1);

    // Check if x_sol_ is empty
    // writeMatrixToFile(x_sol_N_, basePath + "x_sol_N_.txt");

    // Perform matrix operations for all states
    DM ATA = mtimes(transpose(A), A);
    DM ATx = mtimes(transpose(A), ukAll_);  // Transpose x_sol for correct dimension alignment

    // Solve the system to find Bezier coefficients
    bez_GRF = solve(ATA, ATx);

    bez_GRF_ = bez_GRF.T();

    // Output the Bezier Coefficients Matrix
    // writeMatrixToFile(bez_alpha_, basePath + "bez_alpha.txt");
    // std::cout << "Bezier Coefficients Matrix:" << std::endl;
    // std::cout << bez_alpha_ << std::endl;
    // std::cin.get();
}

void NMPC::eigenToCasadi(const Eigen::MatrixXd& vec, casadi::DM& mat) {
    // Resize the CasADi DM matrix to match the dimensions of the Eigen vector
    mat = casadi::DM::zeros(vec.size());

    // Iterate over all elements to copy them
    for (int i = 0; i < vec.size(); ++i) {
        mat(i) = vec(i);  // Assign Eigen value to CasADi DM
    }
}

void NMPC::eigenToCasadi(const Eigen::Vector4d& vec, casadi::DM& mat) {
    // Resize the CasADi DM matrix to match the dimensions of the Eigen vector
    mat = casadi::DM::zeros(vec.size());

    // Iterate over all elements to copy them
    for (int i = 0; i < vec.size(); ++i) {
        mat(i) = vec(i);  // Assign Eigen value to CasADi DM
    }
}


void NMPC::casadiToEigen(const casadi::DM& mat1, Eigen::MatrixXd& mat2) {
    // Resize the Eigen matrix to match the dimensions of the CasADi matrix
    mat2.resize(mat1.rows(), mat1.columns());

    // Iterate over all elements to copy them
    for (int i = 0; i < mat1.rows(); ++i) {
        for (int j = 0; j < mat1.columns(); ++j) {
            mat2(i, j) = static_cast<double>(mat1(i, j));  // Cast CasADi double to Eigen double
        }
    }
}

void NMPC::updateState(double* qin, double* dqin, int* ind, Eigen::Matrix<double, 3, 4> toePos, Eigen::Vector4d state_vec, Eigen::VectorXd GRF_sol_other, const int domain, const int control_tick) { 
    memcpy(q,qin,18*sizeof(double)); 
    memcpy(dq,dqin,18*sizeof(double)); 
    memcpy(contactInd, ind,4*sizeof(int));
    toePos_ = toePos;
    state_other_ = state_vec;

    // std::cout << "GRF_sol_other: " << GRF_sol_other.transpose() << std::endl; std::cin.get();
    GRF_sol_other_.setZero(12, 1); // Initialize to zero
    GRF_sol_other_ = GRF_sol_other;
    
    if(domain > domain_)
    // if(control_tick % 200 == 0)
    {
        isDomainChanged_ = true;
    }
    else
    {
        isDomainChanged_ = false;
    }
    domain_ = domain;
    control_tick_ = control_tick;
    // std::cout << "OTHER AGENT's STATE (" << agent_id_ << "):   " << state_other_.transpose() << std::endl;
} 

void printObstacleMap(const Eigen::MatrixXd& obstacle_map) {
    for (int i = 0; i < obstacle_map.rows(); ++i) {
        for (int j = 0; j < obstacle_map.cols(); ++j) {
            char cell = obstacle_map(i, j) ? 'O' : '-';
            std::cout << cell;
        }
        std::cout << std::endl; // Move to the next line after printing each row
    }
}
//    int ramp_up_iterations = config["parameters"]["ramp_up_iterations"].as<int>();

// if (k == 1)
// {
//     p_g = q.block(4 * 0, i, 2, 1);
// }
void NMPC::generateReferenceTrajectory()
{
    // New way to load parameters from a YAML file, avoids compiling each time a parameter is changed
    double c = config["parameters"]["c"].as<double>();
    double m = config["parameters"]["m"].as<double>();
    double epsilon = config["parameters"]["epsilon"].as<double>();
    double sigma = config["parameters"]["sigma"].as<double>();
    double dth = config["parameters"]["dth"].as<double>();
    double alpha = config["parameters"]["alpha"].as<double>();
    double eta = config["parameters"]["eta"].as<double>();
    double dmin = config["parameters"]["dmin"].as<double>();

    int loopSize = SIM_TIME / params.dt;
    size_t num_obs = Pobs.cols();
    const size_t ramp_up_iterations = 10;

    double Ts = 0.01;

    Eigen::MatrixXd Ae; Ae.setZero(4 * NUMBER_OF_AGENTS, 4 * NUMBER_OF_AGENTS);

    Eigen::MatrixXd A(4, 4);
    Eigen::MatrixXd B(4, 2);
    A << 0, 0, 1, 0,
         0, 0, 0, 1,
         0, 0, -c/m, 0,
         0, 0, 0, -c/m;
    B << 0, 0,
         0, 0,
         1/m, 0,
         0, 1/m;

    // Discretize
    Eigen::MatrixXd Ad = Eigen::MatrixXd::Identity(4, 4) + Ts * A;
    Eigen::MatrixXd Bd = Ts * B;

    for (int agent = 0; agent < NUMBER_OF_AGENTS; ++agent) {
        int blockStartIndex = agent * 4;
        Ae.block(blockStartIndex, blockStartIndex, 4, 4) = Ad;
    }

    Eigen::MatrixXd Be; Be.setZero(4 * NUMBER_OF_AGENTS, 2 * NUMBER_OF_AGENTS);

    for (int agent = 0; agent < NUMBER_OF_AGENTS; ++agent) {
        int blockRowStartIndex = agent * 4;
        int blockColStartIndex = agent * 2;
        Be.block(blockRowStartIndex, blockColStartIndex, 4, 2) = Bd;
    }

    // Define map size and resolution based on operating space
    double x_min = -1.0, x_max = 11.0, y_min = -6.0, y_max = 6.0;
    double resolution = 0.1; // Adjust the resolution for finer grid
    int map_size_x = static_cast<int>((x_max - x_min) / resolution);
    int map_size_y = static_cast<int>((y_max - y_min) / resolution);

    // Initialize the obstacle map
    Eigen::MatrixXd obstacle_map = Eigen::MatrixXd::Zero(map_size_x, map_size_y);

    // Radius around each obstacle to be marked as obstacle
    double obstacle_radius = 0.6;
    int obstacle_radius_cells = static_cast<int>(obstacle_radius / resolution);

    // Map obstacles to the grid
    for (size_t j = 0; j < num_obs; ++j) {
        int obs_x = static_cast<int>((Pobs(0, j) - x_min) / resolution);
        int obs_y = static_cast<int>((Pobs(1, j) - y_min) / resolution);

        if (obs_x >= 0 && obs_x < obstacle_map.rows() && obs_y >= 0 && obs_y < obstacle_map.cols()) {
            for (int dx = -obstacle_radius_cells; dx <= obstacle_radius_cells; ++dx) {
                for (int dy = -obstacle_radius_cells; dy <= obstacle_radius_cells; ++dy) {
                    int nx = obs_x + dx;
                    int ny = obs_y + dy;
                    if (nx >= 0 && nx < obstacle_map.rows() && ny >= 0 && ny < obstacle_map.cols() && (dx * dx + dy * dy <= obstacle_radius_cells * obstacle_radius_cells)) {
                        obstacle_map(nx, ny) = 1; // Mark obstacles on the map
                    }
                }
            }
        }
    }

    // Print obstacle map for debugging
    std::cout << "Obstacle Map:" << std::endl;
    printObstacleMap(obstacle_map);

    // Initialize A* path planning
    std::vector<std::vector<Eigen::Vector2d>> global_paths(NUMBER_OF_AGENTS);
    for (int k = 0; k < NUMBER_OF_AGENTS; ++k) {
        int start_x = static_cast<int>((Pstart_(2 * k) - x_min) / resolution);
        int start_y = static_cast<int>((Pstart_(2 * k + 1) - y_min) / resolution);
        int goal_x = static_cast<int>((GOAL_X - x_min) / resolution);
        int goal_y = static_cast<int>((GOAL_Y - y_min) / resolution);

        std::cout << "Agent " << k << " start: (" << start_x << ", " << start_y << "), goal: (" << goal_x << ", " << goal_y << ")\n";

        std::vector<Eigen::Vector2d> raw_path = a_star(start_x, start_y, goal_x, goal_y, obstacle_map);
        global_paths[k] = filter_waypoints(raw_path);
        if (global_paths[k].empty()) {
            std::cerr << "A* failed to find a path for agent " << k << std::endl;
            return;
        }
        // Convert waypoints back to original resolution
        for (auto& waypoint : global_paths[k]) {
            waypoint[0] = waypoint[0] * resolution + x_min;
            waypoint[1] = waypoint[1] * resolution + y_min;
        }
        // Log the waypoints
        std::cout << "Waypoints for agent " << k << ":\n";
        for (const auto& waypoint : global_paths[k]) {
            std::cout << "(" << waypoint[0] << ", " << waypoint[1] << ")\n";
        }
    }

    // Reference Loop
    Eigen::MatrixXd q; q.setZero(4 * NUMBER_OF_AGENTS, loopSize + 1);

    // Initialize starting positions and velocities
    for (int k = 0; k < NUMBER_OF_AGENTS; ++k) {
        q(4 * k, 0) = Pstart_(2 * k);
        q(4 * k + 1, 0) = Pstart_(2 * k + 1);
        q(4 * k + 2, 0) = 0;
        q(4 * k + 3, 0) = 0;
    }

    // Threshold distance to update waypoint
    double waypoint_threshold = 0.2;

    // Waypoint index for each agent
    std::vector<size_t> waypoint_indices(NUMBER_OF_AGENTS, 0);
    // waypoint_indices[1] = 2;

    for (size_t i = 0; i < loopSize; i++) {
        Eigen::MatrixXd F; F.setZero(2 * NUMBER_OF_AGENTS, 1);

        for (size_t k = 0; k < NUMBER_OF_AGENTS; k++) {
            Eigen::Vector2d qk_p = q.block(4 * k, i, 2, 1);

            // // Get the current waypoint as the goal
            // size_t waypoint_index = std::min(i, global_paths[k].size() - 1);
            // Eigen::Vector2d p_g = global_paths[k][waypoint_index];
            // double d_goal = (qk_p - p_g).norm();

            // Get the current waypoint as the goal
            // Only use the path for 1st agent for coherence (replace with k for multiple agents)
            Eigen::Vector2d p_g = global_paths[0][waypoint_indices[k]];
            double d_goal = (qk_p - p_g).norm();

            // Update the current waypoint based on the threshold distance
            if (d_goal < waypoint_threshold && waypoint_indices[k] < global_paths[0].size() - 1) {
                waypoint_indices[k]++;
                p_g = global_paths[0][waypoint_indices[k]];
                d_goal = (qk_p - p_g).norm();
            }

            // if (k == 1)
            // {
            //     p_g = q.block(4 * 0, i, 2, 1);
            // }

            Eigen::Vector2d F_att; F_att.setZero();
            if (d_goal > 0.01) {
                F_att = -alpha * (qk_p - p_g).normalized();
            } else {
                F_att.setZero();
            }

            Eigen::Vector2d F_rep = Eigen::Vector2d::Zero();
            for (size_t j = 0; j < num_obs; j++) {
                Eigen::Vector2d obstacle(Pobs(0, j), Pobs(1, j));
                double d_obs = (qk_p - obstacle).norm();
                if (d_obs < dmin) {
                    F_rep += eta * (1 / d_obs - 1 / dmin) / (d_obs * d_obs) * ((qk_p - obstacle).normalized());
                }
            }

            Eigen::Vector2d F_agent = Eigen::Vector2d::Zero();
            for (size_t j = 0; j < NUMBER_OF_AGENTS; j++) {
                if (k != j) {
                    Eigen::Vector2d qk_p_other = q.block(4 * j, i, 2, 1);
                    double d_agent = (qk_p - qk_p_other).norm();
                    if (d_agent > 0) {
                        F_agent -= 4 * epsilon * ((6 * pow(sigma, 6) / pow(d_agent, 7)) - (12 * pow(sigma, 12) / pow(d_agent, 13))) * (qk_p - qk_p_other).normalized();
                    }
                }
            }

            F.block(2 * k, 0, 2, 1) = F_att + F_rep + F_agent;
            if (d_goal < 0.001) {
                F.block(2 * k, 0, 2, 1) = Eigen::Vector2d::Zero();
            }

            // Update velocity
            q(4 * k + 2, i + 1) = F(2 * k, 0) / m;
            q(4 * k + 3, i + 1) = F(2 * k + 1, 0) / m;

            // Log positions and forces for debugging
            std::cout << "Agent " << k << ", Step " << i << ": Pos(" << qk_p(0) << ", " << qk_p(1) << "), Force(" << F(0 + k * 2, 0) << ", " << F(1 + k * 2, 0) << ")" << std::endl;
        }

        double scale_factor = (i < ramp_up_iterations) ? (static_cast<double>(i) / ramp_up_iterations) : 1.0;
        q.col(i + 1) = Ae * q.col(i) + Be * F * scale_factor;
    }

    size_t resample_freq = 1;

    for (size_t i = 0; i < loopSize / resample_freq; i++) {
        q.col(i) = q.col(resample_freq * i);
    }

    Eigen::MatrixXd Pr; Pr.setZero(2 * NUMBER_OF_AGENTS, loopSize / resample_freq);
    Eigen::MatrixXd Prd; Prd.setZero(2 * NUMBER_OF_AGENTS, loopSize / resample_freq);

    for (int agent = 0; agent < NUMBER_OF_AGENTS; ++agent) {
        int baseRowPr = agent * 2;
        int baseRowPrd = agent * 2;
        int baseRowQ = agent * 4;

        Pr.block(baseRowPr, 0, 1, loopSize / resample_freq) = q.block(baseRowQ, 1, 1, loopSize / resample_freq);
        Pr.block(baseRowPr + 1, 0, 1, loopSize / resample_freq) = q.block(baseRowQ + 1, 1, 1, loopSize / resample_freq);
        Prd.block(baseRowPrd, 0, 1, loopSize / resample_freq) = q.block(baseRowQ + 2, 1, 1, loopSize / resample_freq);
        Prd.block(baseRowPrd + 1, 0, 1, loopSize / resample_freq) = q.block(baseRowQ + 3, 1, 1, loopSize / resample_freq);
    }


    // // ================================= CUSTOM REFERENCE ========================================//
    // double x0 = 0, y0 = 0, x1 = 0, y1 = -1.0;
    // double scale_factor = 3; // (i < ramp_up_iterations) ? (static_cast<double>(i) / ramp_up_iterations) : 1.0;
    // double time_step = 0.01; // 10ms per iteration
    // double distance_threshold = 1.0; // 1 meter distance for switching y-velocity
    // double vx = 0.15 * scale_factor;
    // double vy0 = -0.15 * scale_factor;
    // double vy1 = 0.15 * scale_factor;

    // for (int i = 0; i < Pr.cols(); i++)
    // {
    //     // Check if we need to switch the y direction for each agent
    //     if (static_cast<int>(x0) / static_cast<int>(distance_threshold) % 2 == 1)
    //     {
    //         vy0 = 0.15 * scale_factor;
    //         vy1 = -0.15 * scale_factor;
    //     }
    //     else
    //     {
    //         vy0 = -0.15 * scale_factor;
    //         vy1 = 0.15 * scale_factor;
    //     }

    //     // Integrate velocities to get positions
    //     x0 += vx * time_step;
    //     y0 += vy0 * time_step;
    //     x1 += vx * time_step;
    //     y1 += vy1 * time_step;

    //     Pr.block(0, i, 4, 1) << x0, 
    //                             y0, 
    //                             x1, 
    //                             y1;

    //     Prd.block(0, i, 4, 1) << vx, 
    //                             vy0, 
    //                             vx, 
    //                             vy1;
    // }



    Pr_refined_ = Pr;
    Prd_refined_ = Prd;

    // Logging Waypoints
    std::filesystem::path baseDir = std::filesystem::current_path().parent_path();
    std::filesystem::path outputsDir = baseDir / "matlab_scripts";

    if (!std::filesystem::exists(outputsDir)) {
        std::filesystem::create_directories(outputsDir);
    }

    std::filesystem::path filePath = outputsDir / "Waypoints.txt";
    std::ofstream waypointFile(filePath);
    if (!waypointFile.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
    } else {
        for (int k = 0; k < NUMBER_OF_AGENTS; ++k) {
            for (const auto& waypoint : global_paths[k]) {
                waypointFile << k << "," << waypoint[0] << "," << waypoint[1] << "\n";
            }
        }
        waypointFile.close();
    }

    // =========================================================================//
    // ================================= LOGGING ===============================//
    // =========================================================================//
    // Step back one directory from the current working directory
    outputsDir = baseDir / "matlab_scripts";

    // Check if the Sim_Outputs directory exists, and create it if it doesn't
    if (!std::filesystem::exists(outputsDir)) {
        std::filesystem::create_directories(outputsDir);
    }

    // Log HL Path
    filePath = outputsDir / "HLPath.txt";
    std::fstream myfile5;
    myfile5.open(filePath, std::fstream::out);
    if (!myfile5.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        // Handle error, possibly exit the function or program
    } else {
        myfile5 << Pr_refined_;
        myfile5.close();
    }

    // Log HL Velocity
    filePath = outputsDir / "HLVelocity.txt";
    std::fstream myfile6;
    myfile6.open(filePath, std::fstream::out);
    if (!myfile6.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        // Handle error, possibly exit the function or program
    } else {
        myfile6 << Prd_refined_;
        myfile6.close();
    }

    // Log Pstart
    filePath = outputsDir / "Pstart.txt";
    std::fstream Pstartf;
    Pstartf.open(filePath, std::fstream::out);
    if (!Pstartf.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        // Handle error, possibly exit the function or program
    } else {
        Pstartf << Pstart_;
        Pstartf.close();
    }

    // Log Pobs
    filePath = outputsDir / "Pobs.txt";
    std::fstream Pobsf;
    Pobsf.open(filePath, std::fstream::out);
    if (!Pobsf.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        // Handle error, possibly exit the function or program
    } else {
        Pobsf << Pobs;
        Pobsf.close();
    }

    // Log Pobs_real
    filePath = outputsDir / "Pobs_real.txt";
    std::fstream Pobs_realf;
    Pobs_realf.open(filePath, std::fstream::out);
    if (!Pobs_realf.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        // Handle error, possibly exit the function or program
    } else {
        Pobs_realf << Pobs_real;
        Pobs_realf.close();
    }
}

void NMPC::setPstart(const Eigen::Matrix<double, 2*NUMBER_OF_AGENTS, 1>& Pstart)
{
    Pstart_ << Pstart;
    agent_Initial_(0) = Pstart_(2*agent_id_);
    agent_Initial_(1) = Pstart_(2*agent_id_+1);

    // params.p_hip(0,0) =  0.183+agent_Initial_(0);  params.p_hip(1,0) = -0.1321+agent_Initial_(1);  params.p_hip(2,0) = 0.01;
    // params.p_hip(0,1) =  0.183+agent_Initial_(0);  params.p_hip(1,1) =  0.1321+agent_Initial_(1);  params.p_hip(2,1) = 0.01;
    // params.p_hip(0,2) = -0.183+agent_Initial_(0);  params.p_hip(1,2) = -0.1321+agent_Initial_(1);  params.p_hip(2,2) = 0.01;
    // params.p_hip(0,3) = -0.183+agent_Initial_(0);  params.p_hip(1,3) =  0.1321+agent_Initial_(1);  params.p_hip(2,3) = 0.01;

    // params.pf(0,0) = toePos_(0, 0) + agent_Initial_(0);  params.pf(1,0) = toePos_(1, 0) + agent_Initial_(0); params.pf(2,0) = toePos_(2, 0);
    // params.pf(0,1) = toePos_(0, 1) + agent_Initial_(0);  params.pf(1,1) = toePos_(1, 1) + agent_Initial_(0); params.pf(2,1) = toePos_(2, 1);
    // params.pf(0,2) = toePos_(0, 2) + agent_Initial_(0);  params.pf(1,2) = toePos_(1, 2) + agent_Initial_(0); params.pf(2,2) = toePos_(2, 2);
    // params.pf(0,3) = toePos_(0, 3) + agent_Initial_(0);  params.pf(1,3) = toePos_(1, 3) + agent_Initial_(0); params.pf(2,3) = toePos_(2, 3);

    // mpc_state_alpha_buffer_<< agent_Initial_(0), 0, agent_Initial_(1), 0;
}

void NMPC::setPobs(Eigen::Matrix<double, 2,NUMBER_OF_OBS>& Pobs_in)
{
    Pobs << Pobs_in;
}

void NMPC::setPobs_real(Eigen::Matrix<double, 2,NUMBER_OF_OBS>& Pobs_in)
{
    Pobs_real << Pobs_in;
}

void NMPC::setAgentID(size_t agent_id)
{
    agent_id_ = agent_id;

    std::string basePath = "tmp/";
    xoutFile_.open(basePath + "xout_" + std::to_string(agent_id_) + ".txt");
    trajFile_.open(basePath + "traj_" + std::to_string(agent_id_) + ".txt");
    uoutFile_.open(basePath + "uout_" + std::to_string(agent_id_) + ".txt");
    feetFile_.open(basePath + "feet_" + std::to_string(agent_id_) + ".txt");
    toutFile_.open(basePath + "tout_" + std::to_string(agent_id_) + ".txt");
    cbfFile_.open(basePath + "cbf_" + std::to_string(agent_id_) + ".txt");
    nmpcSolveTimeFile_.open(basePath + "NMPC_solve_time_" + std::to_string(agent_id_) + ".txt");

}

int NMPC::findNearestIndex() {
    int nearestIndex = 0;
    double minDistance = std::numeric_limits<double>::max();

    Eigen::Vector2d qk; qk << (double)qk_(0), (double)qk_(1);

    for (int i = 0; i < Pr_refined_.cols(); ++i) {
        // Extract the agent's position from Pr_refined_ for the current column
        Eigen::Vector2d agentPosition = Pr_refined_.block(agent_id_ * 2, i, 2, 1);

        // Calculate the Euclidean distance
        double distance = (agentPosition - qk).norm();
        if (distance < minDistance) {
            minDistance = distance;
            nearestIndex = i;
        }
    }

    return nearestIndex;
}


void NMPC::find_nearest_obs() {
    if (Pobs_real.cols() == 0) return; // No obstacles to consider

    // Extract agent position
    double x = static_cast<double>(qk_(0));
    double y = static_cast<double>(qk_(1));

    // Initial settings for finding the nearest obstacle
    double min_distance = std::numeric_limits<double>::max();
    int nearest_index = -1;

    // Loop through all obstacles
    for (int i = 0; i < Pobs_real.cols(); ++i) {
        double obs_x = Pobs_real(0, i);
        double obs_y = Pobs_real(1, i);
        double distance = std::sqrt(std::pow(obs_x - x, 2) + std::pow(obs_y - y, 2));

        if (distance < min_distance) {
            min_distance = distance;
            nearest_index = i;
        }
    }

    // Check the distance to the other agent
    double other_agent_x = static_cast<double>(state_other_(0));
    double other_agent_y = static_cast<double>(state_other_(1));
    double other_agent_distance = std::sqrt(std::pow(other_agent_x - x, 2) + std::pow(other_agent_y - y, 2));

    int due_to_agent = 0;
    if (other_agent_distance < min_distance) {
        min_distance = other_agent_distance;
        obs_near_ << other_agent_x, other_agent_y;
        due_to_agent = 1;
        // cbf_active = 1;
    } else if (nearest_index != -1) {
        obs_near_ = Pobs_real.col(nearest_index);
        due_to_agent = 0;
        // cbf_active = 0;
    }

    if ((min_distance < 0.6 || params.mpc_iter < 150) && due_to_agent == 1)
    {
        cbf_active = 1;
    }
    else
    {
        cbf_active = 0;
    }
    
}

Eigen::Vector4d NMPC::get_lastState() {
    Eigen::Vector4d last_state; 
    last_state << (double)qk_(0), (double)qk_(1), (double)qk_(3), (double)qk_(4); //q[0], q[1], dq[0], dq[1];
    return last_state;
}