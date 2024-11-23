#include "LocoWrapper.hpp"
#include "iostream"
#include "dec_vars_constr_cost.h"
#include <ifopt/ipopt_solver.h>
#include <ifopt/snopt_solver.h>
#include <ifopt/problem.h>

#define PHASE_F 0.15 // in ms - 150 ms of MPC horizon -> N*dt = 15*0.01

using std::cin;
using std::cout;
using namespace ifopt;

void LocoWrapper::stepSetup(double totalStepNum, double step_X, double step_Y, double body_R, double body_P, double body_Y, Eigen::MatrixXd agent_Initial){
    totalStepNum_ = totalStepNum;
    // step_X_ = step_X;
    // step_Y_ = step_Y;
    body_R_ = body_R;
    body_P_ = body_P;
    body_Y_ = body_Y;
    agent_Initial_ = agent_Initial;
}

LocoWrapper::LocoWrapper(int argc, char *argv[]) : Parameters(argc,argv){
    std::string exp_name = "";
    static int counter_agents = 0;
    // std::string filename = "///home/basit/workspace/SRB_NMPC/matlab_scripts/LL/";
    // filename += "agent_" + std::to_string(counter_agents); filename += ".txt"; log_ID = counter_agents;
    // data = std::unique_ptr<DataLog>( new DataLog(filename) ); // make_unique DNE in c++11
    counter_agents++;
    if (counter_agents = 2) { counter_agents = 0; }
    quad = new RobotModel();
    conEst = new ContactEst();
    LL = new LowLevelCtrl();
    VC = new VirtualConstraints();
    PP = new MotionPlanner();

    state = quad->getStatePointer();
    dyn = quad->getDynamicsPointer();
    kin = quad->getKinematicsPointer();
    con = conEst->getConInfoPointer();
    traj = PP->getTrajInfoPointer();
    vcon = VC->getVCPointer();
    ll = LL->getllPointer();
    
    locoTick = 0;
    maxPhase = 1.05;

    onegaitCycle_ = 4; //const
    ts_OptTick_ = TSOPTTICK;
    gridNum_ = NDOMAIN;
    agent_id_ = 0;  // CHANGE FOR MULTIPLE AGENTS
    gaitDomain_ = 0;     // number of current Index column in totalCycleIndex: starting from 0th column =[1,1,1,1]
    body_R_ = BODYROLL;
    body_P_ = BODYPITCH;
    body_Y_ = BODYYAW;
    totalStepNum_ = TOTALSTEPNUM;

    mpc_state_alpha_buffer_.setZero(4,1);

    comTraj_[2] = {0};
    dcomTraj_[2] = {0};
    ddcomTraj_[2] = {0};

    //Eigen::Matrix<double, 2*NUMBER_OF_AGENTS, 1> Pstart;
    //Pstart << -1.5+3,3-3,-1.5+3,3-3+1,-1+3,3-3,-1+3,2-3,-2+3,2-3,0,0,0,0,0,0,0,0,0,0;
    //Pstart << -3,3,-2,3,-1,3,-1,2,-2,2,0,0,0,0,0,0,0,0,0,0;
    //Pstart << 2.0,-2,  2.0,-3;
    domain_ = 0;
    gaitDomain_ = domain_;
}

LocoWrapper::~LocoWrapper(){
    delete quad;
    delete conEst;
    delete LL;
    delete VC;
    delete PP;
}

void LocoWrapper::totalCycleIndexwHalf(size_t gait, size_t gaitCycleNum){

    oneCycleIndex(gait);
    totalCycleIndex_.setOnes(4,4*gaitCycleNum+4);
    for(size_t i=0; i< gaitCycleNum; i++){
        totalCycleIndex_.block(0,i*4+1,4,4)=oneCycleIndex_;
    }
    totalCycleIndex_.block(0,gaitCycleNum*4+1,4,2)=oneCycleIndex_.block(0,0,4,2);
    //totalStepNum_ = totalCycleIndex_.cols();
}

void LocoWrapper::oneCycleIndex(size_t gait){ //total half forward
    if(gait == STAND){
        oneCycleIndex_.setOnes(4,4);
    }    
    if(gait == WALK){
        oneCycleIndex_.setOnes(4,4);
        Eigen::Vector4d leg0stride;
        Eigen::Vector4d leg1stride;
        Eigen::Vector4d leg2stride;
        Eigen::Vector4d leg3stride;
        leg0stride<< 0,1,1,1;
        leg1stride<< 1,0,1,1;
        leg2stride<< 1,1,0,1;
        leg3stride<< 1,1,1,0;
        oneCycleIndex_.block(0,0,4,1) = leg0stride; //half stridekaveh
        oneCycleIndex_.block(0,2,4,1) = leg1stride; //half stride
        oneCycleIndex_.block(0,3,4,1) = leg2stride; //half stride
        oneCycleIndex_.block(0,1,4,1) = leg3stride; //half stride
    }
    if(gait == TROT){
        oneCycleIndex_.setOnes(4,4);
        Eigen::Vector4d leg03stride;
        Eigen::Vector4d leg12stride;
        leg03stride<< 0,1,1,0;
        leg12stride<< 1,0,0,1;
        oneCycleIndex_.block(0,0,4,1) = leg03stride; //half stride
        oneCycleIndex_.block(0,1,4,1) = leg12stride; //full stride
        oneCycleIndex_.block(0,2,4,1) = leg03stride; //full stride
        oneCycleIndex_.block(0,3,4,1) = leg12stride; //half stride       
    }
}

void LocoWrapper::totalCycleIndex(size_t gait, size_t gaitCycleNum){
    oneCycleIndex(gait);
    totalCycleIndex_.setOnes(4, 4*gaitCycleNum+2);
    for(size_t i=0; i< gaitCycleNum; i++){
        totalCycleIndex_.block(0,i*4+1,4,4)=oneCycleIndex_;
    }
    //totalStepNum_ = totalCycleIndex_.cols();
}

void LocoWrapper::plannedCycleIndex(size_t gait){
    size_t cycleNum;
    size_t remainder;
    remainder = (totalStepNum_-2) % onegaitCycle_;
    cycleNum = (totalStepNum_-2-remainder) / onegaitCycle_;

    if (remainder == 0){
        totalCycleIndex(gait, cycleNum);
    }
    else if(remainder == 2){
        totalCycleIndexwHalf(gait, cycleNum);
    }
    else{
        std::cout << "index generation error" << std::endl;
    }
}

void LocoWrapper::calcTau(const double q[18], const double dq[18], const double R[9], const int force[4], size_t gait, size_t ctrlTick, int* checkrunMPC){
    //phaseVar = getPhase(1.0*locoTick, 0.0, 1.0*traj->domLen);   // update phase variable
    // phaseVar = getPhase(1.0*locoTick, 0.0, 1.0*ts_OptTick_*gridNum_);


    if(gait==STAND){
        
        phaseVar = getPhase(1.0*locoTick, 0.0, 1.0*traj->domLen);   // update phase variable
    }else{
        phaseVar = getPhase(1.0*locoTick, 0.0, 1.0*traj->domLen); // 15*0.01  -> N = 15, dt = 0.01 // TODO-INTEGRATION
    }
    quad->updateState(q,dq,R); // update state
    // std::cout << "IS THE DATA ALREADY UPDATED? ===== " << MPC_data_set; std::cout << std::endl;

    static double comx_last = 0;
    static double comy_last = 0;
    
    float footPos[4] = {0};  // DUMMY VARS
    if (gait!=gaitTemp || (phaseVar>maxPhase && gait!=STAND) ){ 
        // Change domain immediately since gait changed
        conEst->forceDomChange();                                       // force con->changeDomain=1 to plan properly
        std::cout << "================================ Time trigger: ================================" << phaseVar << std::endl;
        
        // *checkrunMPC = 1;  MPC_data_available = 0; 
        // Check the footholds
        isDomainChange = 1;
        

        phaseVar = (gaitTemp==STAND) ? 0 : phaseVar;
        getComTrajectory(gait, gaitDomain_);
        // Eigen::Vector4d comTraj; comTraj << comTraj_[0], comTraj_[1], dcomTraj_[0], dcomTraj_[1];
        Eigen::VectorXd comTraj; comTraj.setZero(8);
         comTraj << comTraj_[0], comTraj_[1], comTraj_[2], dcomTraj_[0], dcomTraj_[1], dcomTraj_[2], comOriTraj_[2], dcomOriTraj_[2];
        // Eigen::Vector4d comTraj; comTraj << state->q[0] + 0.001*dcomTraj_[0], state->q[01] + 0.001*dcomTraj_[1], dcomTraj_[0], dcomTraj_[1];
        if (domain_ >= TOTALSTEPNUM) { comTraj(0) = comx_last; comTraj(1) = comy_last; comTraj(3) = 0; comTraj(4) = 0; }
        else {
            comx_last = comTraj(0); comy_last = comTraj(1);
        }
        phaseVar = 0;

        // std::cout << "MPC solution: " << mpc_state_e_x_eventbased_ << std::endl;

        PP->setComDes(comTraj, cbf_active_);
        PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, Pr_refined_, agent_id_, gaitDomain_, mpc_state_e_x_eventbased_);
        conEst->updateConState(footPos,phaseVar,force);
        locoTick = 0;
        gaitDomain_++;
    }else {
        // Wait for impact to change domain
        conEst->updateConState(footPos,phaseVar,force);                 // impact detection
        if (con->changeDomain==1 && gait!=STAND){
            locoTick = 0;
            std::cout << "Contact trigger: " << phaseVar << std::endl;
            phaseVar = 0;
            isDomainChange = 1;
            // *checkrunMPC = 1; MPC_data_available = 0;
            gaitDomain_++;
        }

        getComTrajectory(gait, gaitDomain_);
        if (gaitDomain_ == 0) {
            comTraj_[0] = agent_Initial_(0); comTraj_[1] = agent_Initial_(1);
        }
        cout << "COM (x, y, dx, dy, theta, dtheta) [agent: " << agent_id_ << "]:  " << comTraj_[0] << "  " << comTraj_[1] << "  " << dcomTraj_[0] << "  " << dcomTraj_[1] << "  " << comOriTraj_[2]*(180/3.14) << "  " << dcomOriTraj_[2] << "\n";
        // Eigen::Vector4d comTraj; comTraj << comTraj_[0], comTraj_[1], dcomTraj_[0], dcomTraj_[1];

        Eigen::VectorXd comTraj; comTraj.setZero(8);
        comTraj << comTraj_[0], comTraj_[1], comTraj_[2], dcomTraj_[0], dcomTraj_[1], dcomTraj_[2], comOriTraj_[2], dcomOriTraj_[2];
        // Eigen::Vector4d comTraj; comTraj << state->q[0] + 0.001*dcomTraj_[0], state->q[1] + 0.001*dcomTraj_[1], dcomTraj_[0], dcomTraj_[1];
        // Eigen::Vector4d comTraj; comTraj << 0, 0, 0, 0;
        // comTraj(2) = 0; comTraj(3) = 0;
        if (domain_ >= TOTALSTEPNUM) { comTraj(0) = comx_last; comTraj(1) = comy_last; comTraj(3) = 0; comTraj(4) = 0; }
        else {
            comx_last = comTraj(0); comy_last = comTraj(1);
        }

        if (!STAND) { PP->setComDes(comTraj, cbf_active_);}

        // std::cout << "comTraj: " << comTraj.transpose() << std::endl;

        PP->planTraj(state, kin, conEst, gait, phaseVar, ctrlTick, &motion_params, Pr_refined_, agent_id_, gaitDomain_, mpc_state_e_x_eventbased_);     // plan trajectory
    }

    if (gait != STAND)    // Track the desired GRFs
    {
        getGRFTrajectory(gait, gaitDomain_);
        Eigen::VectorXd GRF_Traj; GRF_Traj.setZero(12);
        for (size_t i = 0; i < 12; i++)
        {
            GRF_Traj(i) = GRF_Traj_[i];
        }
        // VC->setDesiredForce(GRF_Traj);
    }

    quad->updateSwingMatrices(con->ind,con->cnt);                                           // update the jacobians
    VC->updateVirtualConstraints(state, kin, traj, con, gait, phaseVar, &motion_params, ll);     // update VC's
    LL->calcTorque(state, dyn, kin, vcon, con, &ll_params);                             // run low level controller
    data->writeData(state,vcon,traj,ll,ctrlTick,force,phaseVar);                      // log relavent data

    // For SRB NMPC - NMPC runs every 10 ms (100 Hz)  
    dt_ += 0.001;
    if (locoTick % 10 == 0) { 
        dt_ = 0; 
        *checkrunMPC = 1; MPC_data_available = 0;
        phaseVar = 0.0;
    };
    // cout << "locoTick = " << locoTick << "   dt_ = " << dt_ << std::endl;
    locoTick += (ctrlHz)/LL_Hz; // increment locoTick
    gaitTemp = gait; // update the previous gait used
}


void LocoWrapper::logMPC_Data()
{
    std::fstream MPC_Log_Data; ///home/hdsrl7/Desktop/A1_RaiSim_Outputs/exp_
    //MPC_Log_Data.open("/home/hdsrl7/Documents/RAISIM_WORKSPACE/A1_NMPC_oldLL_working/NMPC_Distributed_A1-master/matlab_dbg/MPC_loco.txt",std::fstream::out);
    MPC_Log_Data.open("///home/basit/workspace/SRB_NMPC/matlab_scripts/LL/MPC_loco" + std::to_string(log_ID) + ".txt",std::fstream::out);
    MPC_Log_Data << mpc_state_e_eventbased_;
    MPC_Log_Data.close();

    std::fstream COM_DES;
    COM_DES.open("///home/basit/workspace/SRB_NMPC/matlab_scripts/LL/COM_loco" + std::to_string(log_ID) + ".txt",std::fstream::out);
    COM_DES << com_desired_Traj_eventbased_;
    COM_DES.close();

    std::fstream COP;
    COP.open("///home/basit/workspace/SRB_NMPC/matlab_scripts/LL/COM_loco" + std::to_string(log_ID) + ".txt",std::fstream::out);
    COP << mpc_state_e_u_eventbased_;
    COP.close();
}

void LocoWrapper::setAgentID(size_t agent_id)
{
    agent_id_ = agent_id;

    std::string filename = "///home/basit/workspace/SRB_NMPC/matlab_scripts/LL/";
    filename += "agent_" + std::to_string(agent_id_); filename += ".txt"; log_ID = agent_id_;
    data = std::unique_ptr<DataLog>( new DataLog(filename) ); // make_unique DNE in c++11
}

void LocoWrapper::setPstart(const Eigen::Matrix<double, 2*NUMBER_OF_AGENTS, 1>& Pstart)
{
    Pstart_ << Pstart;
    agent_Initial_(0) = Pstart_(2*agent_id_);
    agent_Initial_(1) = Pstart_(2*agent_id_+1);
    // mpc_state_alpha_buffer_<< agent_Initial_(0), 0, agent_Initial_(1), 0;
}

void LocoWrapper::setPobs(Eigen::Matrix<double, 2,NUMBER_OF_OBS>& Pobs_in)
{
    Pobs << Pobs_in;
}


void LocoWrapper::getComTrajectory(size_t gait, size_t gaitDomain){
    const size_t rownum = 4*TOTALSTEPNUM; //totalstepnumber * 4

    if(gait == STAND){
        comTraj_[0] = 0;
        comTraj_[1] = 0;
        dcomTraj_[0] = 0;
        dcomTraj_[1] = 0;
        ddcomTraj_[0] = 0;
        ddcomTraj_[1] = 0;

        comOriTraj_[0] = 0;
        comOriTraj_[1] = 0;
        comOriTraj_[2] = 0;
        dcomOriTraj_[0] = 0;
        dcomOriTraj_[1] = 0;
        dcomOriTraj_[2] = 0;
        ddcomOriTraj_[0] = 0;
        ddcomOriTraj_[1] = 0;
        ddcomOriTraj_[2] = 0;

        if(gaitDomain >= totalStepNum_-1){
            gaitDomain = totalStepNum_-1;
            Eigen::Matrix<double, 2, 5> alpha_COM_traj_e_current_pos;
            Eigen::Matrix<double, 2, 5> alpha_COM_traj_e_current_vel;

            // SRB Indices: 0 -> x, 1 -> y, 2 -> z, 3 -> xdot, 4 -> ydot, 5 -> zdot

            alpha_COM_traj_e_current_pos.block<1,5>(0,0) = alpha_COM_traj_e_.block(0,0,1,5);
            alpha_COM_traj_e_current_pos.block<1,5>(1,0) = alpha_COM_traj_e_.block(1,0,1,5);
            alpha_COM_traj_e_current_vel.block<1,5>(0,0) = alpha_COM_traj_e_.block(3,0,1,5); // 3 -> xdot
            alpha_COM_traj_e_current_vel.block<1,5>(1,0) = alpha_COM_traj_e_.block(4,0,1,5); // 4 -> ydot

            hzd_dmat* alpha_COM_traj_dmat_pos = new hzd_dmat;
            hzd_dmat* alpha_COM_traj_dmat_vel = new hzd_dmat;
            alpha_COM_traj_dmat_pos->pr = new double[alpha_COM_traj_e_current_pos.rows()*alpha_COM_traj_e_current_pos.cols()];
            alpha_COM_traj_dmat_vel->pr = new double[alpha_COM_traj_e_current_vel.rows()*alpha_COM_traj_e_current_vel.cols()];

            Eigen_mtx_to_hzd_dmat(alpha_COM_traj_e_current_pos, alpha_COM_traj_dmat_pos);
            Eigen_mtx_to_hzd_dmat(alpha_COM_traj_e_current_vel, alpha_COM_traj_dmat_vel);

            HZD_bezier(alpha_COM_traj_dmat_pos, dt_*10, comTraj_); //TODO// just interpolating from the MPC result
            HZD_bezier(alpha_COM_traj_dmat_vel, dt_*10, dcomTraj_); //TODO// little bit weird interpolating. Do it more logically

            delete[] alpha_COM_traj_dmat_pos->pr;
            delete[] alpha_COM_traj_dmat_vel->pr;
            delete alpha_COM_traj_dmat_pos;
            delete alpha_COM_traj_dmat_vel;

            // TODO-Integration of the orientation trajectory
            // comOriTraj_[0] = cubic(dt_*10, 0, PHASE_F,body_R_*(gaitDomain-1),body_R_*(gaitDomain-1),0,0);
            // comOriTraj_[1] = cubic(dt_*10, 0, PHASE_F,body_P_*(gaitDomain-1),body_P_*(gaitDomain-1),0,0);
            // // comOriTraj_[1] = 0;
            // comOriTraj_[2] = cubic(dt_*10, 0, PHASE_F,body_Y_*(gaitDomain-1),body_Y_*(gaitDomain-1),0,0);
            // dcomOriTraj_[0] = cubicDot(dt_*10, 0, PHASE_F,body_R_*(gaitDomain-1),body_R_*(gaitDomain-1),0,0);
            // dcomOriTraj_[1] = cubicDot(dt_*10, 0, PHASE_F,body_P_*(gaitDomain-1),body_P_*(gaitDomain-1),0,0);
            // // dcomOriTraj_[1] = 0;
            // dcomOriTraj_[2] = cubicDot(dt_*10, 0, PHASE_F,body_Y_*(gaitDomain-1),body_Y_*(gaitDomain-1),0,0);
            // ddcomOriTraj_[0] = cubicDotDot(dt_*10, 0, PHASE_F,body_R_*(gaitDomain-1),body_R_*(gaitDomain-1),0,0);
            // ddcomOriTraj_[1] = cubicDotDot(dt_*10, 0, PHASE_F,body_P_*(gaitDomain-1),body_P_*(gaitDomain-1),0,0);
            // ddcomOriTraj_[2] = cubicDotDot(dt_*10, 0, PHASE_F,body_Y_*(gaitDomain-1),body_Y_*(gaitDomain-1),0,0);
        
        }
    }
    
    if(gait != STAND){

        Eigen::Matrix<double, 3, 5> alpha_COM_traj_e_current_pos;
        Eigen::Matrix<double, 3, 5> alpha_COM_traj_e_current_vel;

        Eigen::Matrix<double, 3, 5> alpha_COM_traj_e_current_Ori;
        Eigen::Matrix<double, 3, 5> alpha_COM_traj_e_current_dOri;

        // SRB Indices: 0 -> x, 1 -> y, 2 -> z, 3 -> xdot, 4 -> ydot, 5 -> zdot
        std::cout << "alpha_COM_traj_e: [agent: " << agent_id_ << "]:\n" << alpha_COM_traj_e_ << "\n";

        alpha_COM_traj_e_current_pos.block<1,5>(0,0) = alpha_COM_traj_e_.block(0,0,1,5);
        alpha_COM_traj_e_current_pos.block<1,5>(1,0) = alpha_COM_traj_e_.block(1,0,1,5);
        alpha_COM_traj_e_current_pos.block<1,5>(2,0) = alpha_COM_traj_e_.block(2,0,1,5);

        alpha_COM_traj_e_current_vel.block<1,5>(0,0) = alpha_COM_traj_e_.block(3,0,1,5);
        alpha_COM_traj_e_current_vel.block<1,5>(1,0) = alpha_COM_traj_e_.block(4,0,1,5);
        alpha_COM_traj_e_current_vel.block<1,5>(2,0) = alpha_COM_traj_e_.block(5,0,1,5);

        alpha_COM_traj_e_current_Ori.block<1,5>(0,0) = alpha_COM_traj_e_.block(6,0,1,5);
        alpha_COM_traj_e_current_Ori.block<1,5>(1,0) = alpha_COM_traj_e_.block(7,0,1,5);
        alpha_COM_traj_e_current_Ori.block<1,5>(2,0) = alpha_COM_traj_e_.block(8,0,1,5);
        
        alpha_COM_traj_e_current_dOri.block<1,5>(0,0) = alpha_COM_traj_e_.block(9,0,1,5);
        alpha_COM_traj_e_current_dOri.block<1,5>(1,0) = alpha_COM_traj_e_.block(10,0,1,5);
        alpha_COM_traj_e_current_dOri.block<1,5>(2,0) = alpha_COM_traj_e_.block(11,0,1,5);

        hzd_dmat* alpha_COM_traj_dmat_pos = new hzd_dmat;
        hzd_dmat* alpha_COM_traj_dmat_vel = new hzd_dmat;
        hzd_dmat* alpha_COM_traj_dmat_Ori = new hzd_dmat;
        hzd_dmat* alpha_COM_traj_dmat_dOri = new hzd_dmat;

        alpha_COM_traj_dmat_pos->pr = new double[alpha_COM_traj_e_current_pos.rows()*alpha_COM_traj_e_current_pos.cols()];
        alpha_COM_traj_dmat_vel->pr = new double[alpha_COM_traj_e_current_vel.rows()*alpha_COM_traj_e_current_vel.cols()];

        alpha_COM_traj_dmat_Ori->pr = new double[alpha_COM_traj_e_current_Ori.rows()*alpha_COM_traj_e_current_Ori.cols()];
        alpha_COM_traj_dmat_dOri->pr = new double[alpha_COM_traj_e_current_dOri.rows()*alpha_COM_traj_e_current_dOri.cols()];

        Eigen_mtx_to_hzd_dmat(alpha_COM_traj_e_current_pos, alpha_COM_traj_dmat_pos);
        Eigen_mtx_to_hzd_dmat(alpha_COM_traj_e_current_vel, alpha_COM_traj_dmat_vel);

        Eigen_mtx_to_hzd_dmat(alpha_COM_traj_e_current_Ori, alpha_COM_traj_dmat_Ori);
        Eigen_mtx_to_hzd_dmat(alpha_COM_traj_e_current_dOri, alpha_COM_traj_dmat_dOri);

        if(MPC_data_available){ 
            HZD_bezier(alpha_COM_traj_dmat_pos, dt_*10, comTraj_); //TODO// just interpolating from the MPC result
            HZD_bezier(alpha_COM_traj_dmat_vel, dt_*10, dcomTraj_); //TODO// little bit weird interpolating. Do it more logically
            HZD_bezierd(alpha_COM_traj_dmat_vel, dt_*10, ddcomTraj_); //TODO// little bit weird interpolating. Do it more logically

            HZD_bezier(alpha_COM_traj_dmat_pos, dt_*10, comOriTraj_); //TODO// just interpolating from the MPC result
            HZD_bezier(alpha_COM_traj_dmat_vel, dt_*10, dcomOriTraj_); //TODO// little bit weird interpolating. Do it more logically
            HZD_bezierd(alpha_COM_traj_dmat_vel, dt_*10, ddcomOriTraj_); //TODO// little bit weird interpolating. Do it more logically
        }
        else
        {
            HZD_bezier(alpha_COM_traj_dmat_pos, 0.05, comTraj_); //TODO// just interpolating from the MPC result
            HZD_bezier(alpha_COM_traj_dmat_vel, 0.05, dcomTraj_); //TODO// little bit weird interpolating. Do it more logically
            HZD_bezierd(alpha_COM_traj_dmat_vel, 0.05, ddcomTraj_); //TODO// little bit weird interpolating. Do it more logically

            HZD_bezier(alpha_COM_traj_dmat_pos, 0.05, comOriTraj_); //TODO// just interpolating from the MPC result
            HZD_bezier(alpha_COM_traj_dmat_vel, 0.05, dcomOriTraj_); //TODO// little bit weird interpolating. Do it more logically
            HZD_bezierd(alpha_COM_traj_dmat_vel, 0.05, ddcomOriTraj_); //TODO// little bit weird interpolating. Do it more logically
            // cout << "comx " << comTraj_[0] << "\n"; 
        }

        delete[] alpha_COM_traj_dmat_pos->pr;
        delete[] alpha_COM_traj_dmat_vel->pr;
        delete[] alpha_COM_traj_dmat_Ori->pr;
        delete[] alpha_COM_traj_dmat_dOri->pr;
        delete alpha_COM_traj_dmat_pos;
        delete alpha_COM_traj_dmat_vel;
        delete alpha_COM_traj_dmat_Ori;
        delete alpha_COM_traj_dmat_dOri;

        // comOriTraj_[0] = cubic(dt_*10, 0, PHASE_F,body_R_*(gaitDomain-1),body_R_*gaitDomain,0,0);
        // comOriTraj_[1] = cubic(dt_*10, 0, PHASE_F,body_P_*(gaitDomain-1),body_P_*gaitDomain,0,0);
        // comOriTraj_[2] = cubic(dt_*10, 0, PHASE_F,body_Y_*(gaitDomain-1),body_Y_*gaitDomain,0,0);
        // dcomOriTraj_[0] = cubicDot(dt_*10, 0, PHASE_F,body_R_*(gaitDomain-1),body_R_*gaitDomain,0,0);
        // dcomOriTraj_[1] = cubicDot(dt_*10, 0, PHASE_F,body_P_*(gaitDomain-1),body_P_*gaitDomain,0,0);
        // dcomOriTraj_[2] = cubicDot(dt_*10, 0, PHASE_F,body_Y_*(gaitDomain-1),body_Y_*gaitDomain,0,0);
        // ddcomOriTraj_[0] = cubicDotDot(dt_*10, 0, PHASE_F,body_R_*(gaitDomain-1),body_R_*gaitDomain,0,0);
        // ddcomOriTraj_[1] = cubicDotDot(dt_*10, 0, PHASE_F,body_P_*(gaitDomain-1),body_P_*gaitDomain,0,0);
        // ddcomOriTraj_[2] = cubicDotDot(dt_*10, 0, PHASE_F,body_Y_*(gaitDomain-1),body_Y_*gaitDomain,0,0);

        // cout << "GET COM DBG 12\n";
    }
}

void LocoWrapper::getGRFTrajectory(size_t gait, size_t gaitDomain){
   
    if(gait != STAND){

        Eigen::Matrix<double, 12, 5> alpha_GRF_traj_e_current_pos;
        Eigen::Matrix<double, 12, 5> alpha_GRF_traj_e_current_vel;

        // SRB Indices: 0 -> x, 1 -> y, 2 -> z, 3 -> xdot, 4 -> ydot, 5 -> zdot

        for (size_t i = 0; i < 12; i++)
        {
            alpha_GRF_traj_e_current_pos.block<1,5>(i,0) = alpha_GRF_traj_e_.block(i,0,1,5);
        }
        
        // alpha_GRF_traj_e_current_pos.block<1,5>(1,0) = alpha_GRF_traj_e_.block(1,0,1,5);
        // alpha_GRF_traj_e_current_pos.block<1,5>(2,0) = alpha_GRF_traj_e_.block(2,0,1,5);

        hzd_dmat* alpha_GRF_traj_dmat_pos = new hzd_dmat;
        alpha_GRF_traj_dmat_pos->pr = new double[alpha_GRF_traj_e_current_pos.rows()*alpha_GRF_traj_e_current_pos.cols()];

        Eigen_mtx_to_hzd_dmat(alpha_GRF_traj_e_current_pos, alpha_GRF_traj_dmat_pos);

        if(MPC_data_available){ 
            HZD_bezier(alpha_GRF_traj_dmat_pos, dt_*10, GRF_Traj_); //TODO// just interpolating from the MPC result
        }
        else
        {
            HZD_bezier(alpha_GRF_traj_dmat_pos, 0.05, GRF_Traj_); //TODO// just interpolating from the MPC result
        }

        delete[] alpha_GRF_traj_dmat_pos->pr;
        delete alpha_GRF_traj_dmat_pos;
    }

    // std::cout << "GRFS (z): ";
    // for (size_t i = 1; i <= 4; i++)
    // {
    //     std::cout << GRF_Traj_[i * 3 - 1] << "   ";
    // }
    std::cout << std::endl;
    
}


void LocoWrapper::printSize(const Eigen::MatrixXd& mat)
{
    cout << "\nRows(): " << mat.rows() << "   cols(): " << mat.cols() << std::endl;
    cout << mat << std::endl;
} 

void LocoWrapper::set_MPC_DATA(Eigen::MatrixXd alpha_COM, Eigen::MatrixXd alpha_GRF, Eigen::MatrixXd MPC_sol, bool avail, int cbf_status){
    alpha_COM_traj_e_.setZero(12, 5);
    alpha_GRF_traj_e_.setZero(12, 5);

    alpha_COM_traj_e_ = alpha_COM;
    alpha_GRF_traj_e_ = alpha_GRF;

    mpc_state_e_x_eventbased_ = MPC_sol;
    MPC_data_available = avail;

    cbf_active_ = cbf_status;
};

void LocoWrapper::reset_MPC_DATA()
{
    // DO I NEED TO DO ANYTHING WITH THE alpha_COM_traj_e_ ? or HOW TO ENFORCE LAST known COM positions before MPC solution becomes available? 
}