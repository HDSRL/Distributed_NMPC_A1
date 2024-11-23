#ifndef SHARED_DATA
#define SHARED_DATA

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Sparse"

#include "mutex"
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <vector>
#include <algorithm>

#define SET_DATA 1
#define GET_DATA 0
#define HL_DATA 1
#define LL_DATA 0
#define HORIZON_SIM 12
#define NUMBER_OF_AGENTS 2 // Define the number of agents as per your requirement

boost::mutex mtx;

struct sharedData {
    int resetRun = 0;

    // LL sets, used by HL
    double q[18] = {0};
    double dq[18] = {0};
    int ind[4] = {1};
    int domain = 0;
    int control_tick = 0; 
    int domain_change = 0;
    int cbf_active = 0;

    Eigen::Matrix<double, 3, 4> toePos = Eigen::MatrixXd::Zero(3, 4);

    // HL sets, used by LL
    Eigen::MatrixXd alpha_COM = Eigen::MatrixXd::Zero(12, 5);
    Eigen::MatrixXd alpha_GRF = Eigen::MatrixXd::Zero(12, 5);
    Eigen::MatrixXd MPC_sol_;

    // Both set, both use
    std::vector<Eigen::Vector4d> last_state = std::vector<Eigen::Vector4d>(2, Eigen::Vector4d::Zero()); // Changed to a vector of size 2
    Eigen::MatrixXd GRFs_sol_ = Eigen::MatrixXd::Zero(2, 12); // Changed to a vector of size 2

    int runMPC = 1;
    int MPC_data_available = 0;
};

// Vector of sharedData for multiple agents
std::vector<sharedData> agentData(NUMBER_OF_AGENTS);

void updateData(int setget, int highlow, sharedData* newData, int agent_number) {
    boost::lock_guard<boost::mutex> guard(mtx);

    if (agent_number < 0 || agent_number >= agentData.size()) return;

    sharedData& data = agentData[agent_number];

    if (setget == SET_DATA) {
        // SET DATA
        if (highlow == HL_DATA) {
            data.alpha_COM = newData->alpha_COM;
            data.alpha_GRF = newData->alpha_GRF;
            data.MPC_sol_ = newData->MPC_sol_;
            data.domain = newData->domain;
            data.domain_change = newData->domain_change;
            data.control_tick = newData->control_tick;
            data.MPC_data_available = newData->MPC_data_available;
            data.cbf_active = newData->cbf_active;
            data.runMPC = newData->runMPC;
            data.resetRun = newData->resetRun;
            data.last_state[agent_number] = newData->last_state[agent_number];
            data.GRFs_sol_.row(agent_number) = newData->GRFs_sol_.row(agent_number);
        } else {
            std::copy(newData->q, newData->q + 18, data.q);
            std::copy(newData->dq, newData->dq + 18, data.dq);
            std::copy(newData->ind, newData->ind + 4, data.ind);
            data.domain = newData->domain;
            data.domain_change = newData->domain_change;
            data.control_tick = newData->control_tick;
            data.toePos = newData->toePos;
            if (data.resetRun == -1 && newData->resetRun == -2) {
                data.resetRun = 0;
            }
            if (data.resetRun == 0) {
                data.runMPC = newData->runMPC;
                data.MPC_data_available = newData->MPC_data_available;
                data.cbf_active = newData->cbf_active;
            }
        }
    } else { // GET DATA
        if (highlow == HL_DATA) {
            std::copy(data.q, data.q + 18, newData->q); 
            std::copy(data.dq, data.dq + 18, newData->dq);
            std::copy(data.ind, data.ind + 4, newData->ind);
            newData->MPC_data_available = data.MPC_data_available;
            newData->cbf_active = data.cbf_active;
            newData->runMPC = data.runMPC;
            newData->resetRun = data.resetRun;
            newData->toePos = data.toePos;
            newData->domain = data.domain;
            newData->domain_change = data.domain_change;
            newData->control_tick = data.control_tick;
            newData->last_state[agent_number] = data.last_state[agent_number];
            newData->GRFs_sol_.row(agent_number) = data.GRFs_sol_.row(agent_number);
        } else {
            newData->alpha_COM = data.alpha_COM;
            newData->alpha_GRF = data.alpha_GRF;
            newData->MPC_sol_ = data.MPC_sol_;
            // newData->GRFs_sol_[agent_number] = data.GRFs_sol_[agent_number];
            newData->MPC_data_available = data.MPC_data_available;
            newData->cbf_active = data.cbf_active;
            newData->runMPC = data.runMPC;
            newData->resetRun = data.resetRun;
            newData->domain = data.domain;
            newData->domain_change = data.domain_change;
            newData->control_tick = data.control_tick;
        }
    }
}

void backupData(sharedData& copy_to, const sharedData& copy_from) {
    // Use a lock guard to ensure thread safety if necessary
    boost::lock_guard<boost::mutex> guard(mtx);

    // Copying all fields from copy_from to copy_to
    copy_to.resetRun = copy_from.resetRun;
    
    // For arrays, use std::copy
    std::copy(std::begin(copy_from.q), std::end(copy_from.q), std::begin(copy_to.q));
    std::copy(std::begin(copy_from.dq), std::end(copy_from.dq), std::begin(copy_to.dq));
    std::copy(std::begin(copy_from.ind), std::end(copy_from.ind), std::begin(copy_to.ind));
    
    copy_to.domain = copy_from.domain;
    copy_to.toePos = copy_from.toePos;

    // For Eigen matrices, you can simply use assignment as they support deep copy by default
    copy_to.alpha_COM = copy_from.alpha_COM;
    copy_to.alpha_GRF = copy_from.alpha_GRF;
    copy_to.MPC_sol_ = copy_from.MPC_sol_;
    copy_to.GRFs_sol_.row(0) = copy_from.GRFs_sol_.row(0);
    copy_to.GRFs_sol_.row(1) = copy_from.GRFs_sol_.row(1);
    copy_to.last_state[0] = copy_from.last_state[0];
    copy_to.last_state[1] = copy_from.last_state[1];

    copy_to.runMPC = copy_from.runMPC;
    copy_to.MPC_data_available = copy_from.MPC_data_available;
    copy_to.cbf_active = copy_from.cbf_active;
}

#endif
