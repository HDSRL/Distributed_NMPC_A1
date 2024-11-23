#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include "nmpc_params.h"
#include <array>
#include <vector>

using namespace casadi;
namespace fs = std::filesystem;

int main(){

    // File path setup
    bool code_gen = 1; // Set to 1 to generate code
    std::string basePath = "tmp/";
    fs::create_directories(basePath); // Ensure directory exists

    // =========================== SRB NMPC =========================== //
    Params params(0.01); // 0.1 is the sample time for MPC 
    const double sim_time = 10.0; 
    // const size_t N = 12; // Prediction horizon
    const double dt = 0.01; // Sample time
    // const double Fs = 1.0/dt; // Control Frequency
    // const double mu = 0.6; // Friction Coefficient
    const float sim_length = sim_time/dt; // Simulation length
    // const double height = 0.28; // Stand height
    const int n_states = 12; // Number of states
    const int n_inputs = 12; // Number of inputs
    const int n_controls = 12; // Number of inputs
    const int gait = 8; // Gait type
    const double G = 9.81; // Gravity
    const double mass = 12.4530; // Mass
    // const int agentIndex = 1;

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

    // Define symbolic
    casadi::SXVector footholds_vec = {
        SX::sym("pf1x"), SX::sym("pf1y"), SX::sym("pf1z"),
        SX::sym("pf2x"), SX::sym("pf2y"), SX::sym("pf2z"),
        SX::sym("pf3x"), SX::sym("pf3y"), SX::sym("pf3z"),
        SX::sym("pf4x"), SX::sym("pf4y"), SX::sym("pf4z")
    };

    casadi::SX states = casadi::SX::vertcat(states_vec);
    casadi::SX controls = casadi::SX::vertcat(controls_vec);
    casadi::SX footholds = casadi::SX::vertcat(footholds_vec);

    SX U = SX::sym("U", n_controls, params.N);

    SX P = SX::sym("P", n_states + n_states * params.N + 2 + 16 + 12); // +2 is for nearest obstacle coordinates
                                                                       // +16 is for other agent's (pos) + (vel) + GRFs (order?) : 2 + 2 + 12
                                                                       // +12 for footholds
    SX X = SX::sym("X", n_states, params.N + 1);

    SX slack_var = SX::sym("slack_var", 1, 1);
    SX OPT_variables = vertcat(reshape(X, n_states * (params.N + 1), 1), reshape(U, n_inputs * params.N, 1), slack_var);

    DM qk = DM::zeros(n_states);
    qk(2) = 0.08; // Assuming z position is constant or predefined

    double t0 = 0.0;
    DM xx(n_states, sim_length);
    xx(Slice(0, n_states), 0) = qk;  

    DM u0 = DM::zeros(params.N, params.n_inputs); // Control inputs for N horizon steps, not sim_steps
    
    // States decision variables initialization
    DM X0 = repmat(qk, 1, params.N+1);
    DM xd = repmat(qk, params.N, 1);

    // Write initial symbolic states and controls to file
    // writeMatrixToFile(states, basePath + "states.txt");
    // writeMatrixToFile(controls, basePath + "controls.txt");
    // writeMatrixToFile(qk, basePath + "qk.txt");
    // writeMatrixToFile(xd, basePath + "xd.txt"); 
    // writeMatrixToFile(xx, basePath + "xx.txt");
    // writeMatrixToFile(X, basePath + "X.txt");
    // writeMatrixToFile(P, basePath + "P.txt");
    // writeMatrixToFile(U, basePath + "U.txt");
    // writeMatrixToFile(OPT_variables, basePath + "OPT_variables.txt");
    // writeMatrixToFile(params, basePath + "params.txt");

    double current_time = 0.0;
    size_t mpciter = 0;

    Args args;
    args.lbg = DM::zeros(n_states * (params.N + 1) + 4*4*params.N + params.N, 1);
    args.ubg = DM::zeros(n_states * (params.N + 1) + 4*4*params.N + params.N, 1);

    DM uout = DM::zeros(n_controls, sim_length); // Vector to store control outputs
    DM xout = DM::zeros(n_states, sim_length);
    DM feet = DM::zeros(12, sim_length); // Vector to store foot positions
    DM traj = DM::zeros(12, sim_length); // Vector to store desired trajectory
    DM tout = DM::zeros(1, sim_length);

    // while(mpciter < sim_time / dt) {

        // DEBUG SESSION //
        params.current_time = current_time;
        params.mpc_iter = mpciter;

        
        
        // Updates the xd to obtain desired 12-state reference for next N domains
        planner(current_time, qk, xd, gait, params); // Assuming gait and other details are handled inside

        // writeMatrixToFile(states, basePath + "states.txt");
        // writeMatrixToFile(controls, basePath + "controls.txt");
        // writeMatrixToFile(qk, basePath + "qk.txt");
        // writeMatrixToFile(xd, basePath + "xd.txt"); 
        // writeMatrixToFile(xx, basePath + "xx.txt");
        // writeMatrixToFile(X, basePath + "X.txt");
        // writeMatrixToFile(P, basePath + "P.txt");
        // writeMatrixToFile(U, basePath + "U.txt");
        // writeMatrixToFile(OPT_variables, basePath + "OPT_variables.txt");
        // writeMatrixToFile(params, basePath + "params.txt");

        DM obs = DM::zeros(2, 1); obs(0) = 0.5; obs(1) = 0.1; DM other_agent_info = DM::zeros(16, 1);
        args.p = vertcat(qk, xd, obs, other_agent_info); // test

        // writeMatrixToFile(args.p, basePath + "args_p.txt");
        
        SX t_rot = casadi::SX::vertcat({states(6), states(7), states(8)});  // Extract phi, theta, and psi into a column vector for Rotz
        SX Rz = Rotz(t_rot);    // Call Rotz with the column vector
        SX z3 = casadi::SX::zeros(3, 3);    // Create a 3x3 zero matrix
        SX i3 = casadi::SX::eye(3);         // Create a 3x3 identity matrix

        DM J = DM(params.J); // Assuming a paramse-5; J(2,2) = 0.064713601;

        SX omega = casadi::SX::vertcat({states(9), states(10), states(11)}); // Extracting angular velocities (dphi, dtheta, dpsi)
        SX skew_omega = skewSym(omega);
        SX J_sx = SX(J); // Convert J from DM to SX
        // SX IwI = mtimes(mtimes(mtimes(Rz, J_sx), Rz.T()), skew_omega) * mtimes(mtimes(Rz, J_sx), Rz.T());
        SX IwI = mtimes(mtimes(mtimes(Rz, solve(J_sx, Rz.T())), skew_omega), mtimes(mtimes(Rz, J_sx), Rz.T()));

        // writeMatrixToFile(Rz, basePath + "Rz.txt");
        // writeMatrixToFile(-IwI, basePath + "IwI.txt");
        // writeMatrixToFile(J_sx, basePath + "J_sx.txt");

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
        A(Slice(9, 12), Slice(9, 12)) = -IwI; //-IwI

        // writeMatrixToFile(A, basePath + "A.txt");
        
        SX D = SX::zeros(12,1); // D is a 12x1 vector
        D(Slice(3,6)) = SX::vertcat({SX::zeros(2,1), -G}); // Fill the 4th to 6th rows

        // writeMatrixToFile(D, basePath + "D.txt");

        // SX p_foot_sx = SX(params.pf);
        SX p_foot_sx = reshape(footholds.T(), 3, 4);
        writeMatrixToFile(p_foot_sx, basePath + "p_foot_sx.txt");

        SX com_vector = SX::vertcat({states(0), states(1), states(2)});
        SX com_vector_expanded = repmat(com_vector, 1, 4); // Replicate across columns to match p_foot_sx dimensions

        SX rd = p_foot_sx - com_vector_expanded;
        
        SX skew_rd = skewSym(rd);
        SX B_tmp = mtimes(mtimes(Rz, J_sx), Rz.T()) * skew_rd;
        // writeMatrixToFile(B_tmp, basePath + "B_tmp.txt");
        // writeMatrixToFile(skew_rd, basePath + "skew_rd.txt");
        // writeMatrixToFile(rd, basePath + "rd.txt");

        SX B = SX::vertcat({
            SX::zeros(3, 12),
            SX::repmat(1/mass * i3, 1, 4),
            SX::zeros(3, 12),
            mtimes(mtimes(Rz, J_sx), Rz.T()) * skew_rd}); //mtimes(mtimes(Rz, J_sx), Rz.T()) * skew_rd

        // writeMatrixToFile(B, basePath + "B.txt");

        SX rhs = mtimes(A, states) + mtimes(B, controls) + D;



        // writeMatrixToFile(rhs, basePath + "rhs.txt");

        Function f_dyn = Function("f_dyn", {states, controls, footholds}, {rhs}, {"st", "con", "pf"}, {"rhs"});
        // Declare variables to hold the outputs
        SX g_force, g_ub, g_lb;

        // Call the function
        force_inequality(U, params, g_force, g_ub, g_lb);

        // writeMatrixToFile(g_force, basePath + "g_force.txt");
        // writeMatrixToFile(g_ub, basePath + "g_ub.txt");
        // writeMatrixToFile(g_lb, basePath + "g_lb.txt");
        // writeMatrixToFile(params, basePath + "params.txt");

        SX obj = 0, g_eq; 
        equalityconstraints(X, U, P, params, obj, g_eq, slack_var, f_dyn);

        // writeMatrixToFile(obj, basePath + "obj.txt");
        // writeMatrixToFile(g_eq, basePath + "g_eq.txt");

        

        // args.lbg = DM::zeros(n_states * (N + 1) + 4*4*N, 1);
        // args.ubg = DM::zeros(n_states * (N + 1) + 4*4*N, 1);
        set_statebounds(mpciter, params, args);
        set_constraint_bounds(args, g_lb, g_ub, n_states, params.N);


        SX g_obs;
        // collision_constraints(X, P, params, g_obs, args, slack_var); // g_obs, o_lb, o_ub are being set inside the function 
        cbf_constraints(X, P, params, g_obs, args, slack_var); // g_obs, o_lb, o_ub are being set inside the function

        SX g = vertcat(g_eq, g_force, g_obs);

        // writeMatrixToFile(args.lbx, basePath + "args_lbx.txt");
        // writeMatrixToFile(args.ubx, basePath + "args_ubx.txt");
        writeMatrixToFile(g, basePath + "g.txt");
        writeMatrixToFile(args.lbg, basePath + "args_lbg.txt");
        writeMatrixToFile(args.ubg, basePath + "args_ubg.txt");
        
        // Dictionary to hold NLP problem definitions
        SXDict nlp_prob = {{"f", obj}, {"x", OPT_variables}, {"g", g}, {"p", P}};

        // Define options for the solver
        Dict opts;

        // For ipopt
        // Set Ipopt-specific options
        opts["ipopt.max_iter"] = 10;  // Replace Max_mpciter with its actual value
        opts["ipopt.print_level"] = 0;  // Can be changed to 0 or 3 based on verbosity required
        opts["print_time"] = 0;  // Disable printing solver time
        opts["ipopt.acceptable_tol"] = 1e-3;  // Tolerance for stopping criterion
        opts["ipopt.acceptable_obj_change_tol"] = 1e-3;  // Objective change tolerance for stopping

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
        DM other_agent_IC = DM::zeros(16, 1);
        DM footholds_pf = DM::zeros(12, 1);
        DM x0u0= vertcat(X0_col, u0_col, DM(0));
        args.x0 = x0u0;
        args.p = vertcat(qk, xd, obs, other_agent_IC, footholds_pf); // test
        // Solve the problem
        DMDict arg = {{"x0", args.x0}, {"lbx", args.lbx}, {"ubx", args.ubx}, {"lbg", args.lbg}, {"ubg", args.ubg}, {
            "p", args.p}};
        DMDict res = solver(arg);

        // writeMatrixToFile(res.at("x"), basePath + "solution.txt");
        // writeMatrixToFile(args.p, basePath + "args_p.txt");


        DM sol = DM::zeros(n_states * (params.N + 1) + n_controls * params.N, 1); // sol.x should be initialized with actual solution data
        sol = res.at("x");


        // =========================== CODE GEN =========================== //
        // // Generate C code for the NLP functions
        if (code_gen == 1)
        {
            // // file name
            std::string file_name = "nlp_code";
            // code predix
            std::string prefix_code = fs::current_path().string() + "/";

            casadi::Dict opts = casadi::Dict();
            opts["cpp"] = false;
            opts["with_header"] = false;


            solver.generate_dependencies(file_name + ".c");
            // f.generate_dependencies(file_name + "_f" + ".c");

            CodeGenerator codegen = CodeGenerator(file_name + "_f" + ".c", opts);
            codegen.add(f_dyn);
            // codegen.add(solver);
            // codegen.with_header = true;
            // codegen.verbose = true;
            codegen.generate();

            // shared library prefix
            std::string prefix_lib = fs::current_path().string() + "/";

            // compile c code to a shared library
            // std::string compile_command_solver = "gcc -fPIC -shared -O3 " + 
            //     prefix_code + file_name + ".c -o " +
            //     prefix_lib + file_name + ".so";

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

            code_gen = 0; // do not generate code again
        }
        // if (code_gen == 1)
        // {
        //     // // file name
        //     std::string file_name = "nlp_f";
        //     // code predix
        //     std::string prefix_code = fs::current_path().string() + "/";


        //     // solver.generate_dependencies(file_name + "_solver" + ".c");
        //     // f.generate_dependencies(file_name + "_f" + ".c");

        //     CodeGenerator codegen(file_name + ".c");
        //     codegen.add(f);
        //     // codegen.add(solver);
        //     codegen.generate();

        //     // shared library prefix
        //     std::string prefix_lib = fs::current_path().string() + "/";

        //     // compile c code to a shared library
        //     std::string compile_command_solver = "gcc -fPIC -shared -O3 " + 
        //         prefix_code + file_name + ".c -o " +
        //         prefix_lib + file_name + ".so";

        //     // std::string compile_command_f = "gcc -fPIC -shared -O3 " + 
        //     //     prefix_code + file_name + "_f" + ".c -o " +
        //     //     prefix_lib + file_name + ".so";

        //     std::cout << compile_command_solver << std::endl;
        //     // std::cout << compile_command_f << std::endl;

        //     int compile_flag_solver = std::system(compile_command_solver.c_str());
        //     casadi_assert(compile_flag_solver==0, "Compilation failed");
        //     std::cout << "Compilation successed!" << std::endl;

        //     // int compile_flag_f = std::system(compile_command_f.c_str());
        //     // casadi_assert(compile_flag_f==0, "Compilation failed");
        //     // std::cout << "Compilation successed!" << std::endl;

        //     code_gen = 0; // do not generate code again
        // }
        // =========================== CODE GEN =========================== //



        // Explicitly cast start and stop indices to casadi_int
        casadi_int start_index_u = static_cast<casadi_int>(n_states * (params.N + 1));
        casadi_int end_index = static_cast<casadi_int>(sol.numel()-1);

        // Get controls from the solution
        DM u = reshape(sol(Slice(start_index_u, end_index)).T(), n_controls, params.N).T();

        // Get solution trajectory
        DM sol_x_N = reshape(sol(Slice(0, static_cast<casadi_int>(n_states * (params.N + 1)))).T(), n_states, params.N + 1).T();
        // writeMatrixToFile(sol_x_N, basePath + "x_sol_N.txt");

        // Update logs
        uout(Slice(), mpciter) = u(0, Slice());     // Store current control outputs

        // writeMatrixToFile(u(0, Slice()).T(), basePath + "u_sol.txt");
        // writeMatrixToFile(u.T(), basePath + "u_sol_N.txt");

        // Call the shift function
        // writeMatrixToFile(qk, basePath + "x_sol.txt");
        shift(t0, qk, u0, dt, u, f_dyn);

        // Retrieve solution trajectory (assuming sol.x is available)
        // X0 = reshape(X0.T(), n_states, N + 1).T();
        //X0 = vertcat(sol_x_N(Slice(1, sol_x_N.size1()), Slice()), sol_x_N(Slice(X0.size1() - 1, sol_x_N.size1()), Slice()));
        X0 = vertcat(sol_x_N(Slice(1, sol_x_N.size1()), Slice()), sol_x_N(Slice(sol_x_N.size1() - 1, sol_x_N.size1()), Slice()));
        // writeMatrixToFile(X0, basePath + "x_init_next_step.txt");

        // Data Logging

        xout(Slice(), mpciter) = qk.T(); // Store current states, transposed
        feet(Slice(), mpciter) = reshape(params.pf, 1, 12); // Store reshaped foot positions
        tout(0, mpciter) = t0;          // Store current time

        
        // Function f = Function("f", {states, controls}, {rhs}, {"states", "controls"}, {"rhs"});

        mpciter += 1;
        current_time += dt;
        std::cout << "Current Time: " << current_time << std::endl;
        
        // SX g_force; DM g_ub, g_lb;
        // force_inequality(U, params, g_force, g_ub, g_lb);

    // }

    writeMatrixToFile(xout.T(), basePath + "xout.txt");
    writeMatrixToFile(feet.T(), basePath + "feet.txt");
    writeMatrixToFile(tout.T(), basePath + "tout.txt");
    writeMatrixToFile(uout.T(), basePath + "uout.txt");
    // =========================== EXAMPLE CODE =========================== //

    /*  Test problem
     *
     *    min x0^2 + x1^2
     *    s.t.    x0 + x1 - 10 = 0
     */

    // // Optimization variables
    // casadi::MX x = casadi::MX::sym("x", 2);

    // // Objective
    // casadi::MX f = x(0)*x(0) + x(1)*x(1);

    // // Constraints
    // casadi::MX g = x(0)+x(1)-10;

    // // Create an NLP solver instance
    // casadi::Function solver = casadi::nlpsol("solver", "ipopt", {{"x", x}, {"f", f}, {"g", g}});

    // // file name
    // std::string file_name = "nlp_code";
    // // code predix
    // std::string prefix_code = fs::current_path().string() + "/";

    // // Generate C code for the NLP functions
    // solver.generate_dependencies(file_name + ".c");

    // // shared library prefix
    // std::string prefix_lib = fs::current_path().string() + "/";

    // // compile c code to a shared library
    // std::string compile_command = "gcc -fPIC -shared -O3 " + 
    //     prefix_code + file_name + ".c -o " +
    //     prefix_lib + file_name + ".so";

    // std::cout << compile_command << std::endl;
    // int compile_flag = std::system(compile_command.c_str());
    // casadi_assert(compile_flag==0, "Compilation failed");
    // std::cout << "Compilation successed!" << std::endl;

    // return 0;
}


// int main(){
//     /*  Test problem
//      *
//      *    min x0^2 + x1^2
//      *    s.t.    x0 + x1 - 10 = 0
//      */

//     // Optimization variables
//     casadi::SX x = casadi::SX::sym("x", 2);

//     // Objective
//     casadi::SX f = x(0)*x(0) + x(1)*x(1);

//     // Constraints
//     casadi::SX g = x(0)+x(1)-10;

//     // Create an NLP solver instance
//     casadi::Function solver = casadi::nlpsol("solver", "ipopt", {{"x", x}, {"f", f}, {"g", g}});

//     // file name
//     std::string file_name = "nlp_code";
//     // code predix
//     std::string prefix_code = fs::current_path().string() + "/";

//     // Generate C code for the NLP functions
//     solver.generate_dependencies(file_name + ".c");

//     // shared library prefix
//     std::string prefix_lib = fs::current_path().string() + "/";

//     // compile c code to a shared library
//     std::string compile_command = "gcc -fPIC -shared -O3 " + 
//         prefix_code + file_name + ".c -o " +
//         prefix_lib + file_name + ".so";

//     std::cout << compile_command << std::endl;
//     int compile_flag = std::system(compile_command.c_str());
//     casadi_assert(compile_flag==0, "Compilation failed");
//     std::cout << "Compilation successed!" << std::endl;

//     return 0;
// }