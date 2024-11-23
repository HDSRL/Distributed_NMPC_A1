//
// Authror: Basit M. Imran.
// Date : 2024-04-17
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include <yaml-cpp/yaml.h>
#include "raisim/OgreVis.hpp"
#include "randyImguiPanel.hpp"
#include "raisimKeyboardCallback.hpp"
#include "helper.hpp"
#include "helper2.hpp"
// #include "MPC_dist.hpp"
#include "nmpc.hpp"
#include "OtherUtils.hpp"
#include "timer.h"
#include <cstdlib> // For rand() and srand()
#include <ctime>   // For time()
#include <fstream>
#include <iostream>
//HDSRL header
#include "LocoWrapper.hpp"
#include "shared_structs.hpp"

bool rough_terrain_en = 0;
bool obstacle_en = 1;

using std::cout;
using std::cin;

std::vector<int> first_time(NUMBER_OF_AGENTS, 1);
std::vector<int> sec_time(NUMBER_OF_AGENTS, 1);
std::vector<int> runMPC(NUMBER_OF_AGENTS, 0);
std::vector<sharedData> HLData, LLData, HLData_Backup, LLData_Backup;

void setupCallback() {
    raisim::OgreVis *vis = raisim::OgreVis::get();

    /// light
    vis->getLight()->setDiffuseColour(1, 1, 1);
    vis->getLight()->setCastShadows(false);
    Ogre::Vector3 lightdir(-3,3,-0.5); // Light shines on ROBOTS top/front/right side
    // Ogre::Vector3 lightdir(-3,-3,-0.5); // Light shines on ROBOTS top/front/left side
    lightdir.normalise();
    vis->getLightNode()->setDirection({lightdir});
    vis->setCameraSpeed(300);

    vis->addResourceDirectory(raisim::loadResource("material"));
    vis->loadMaterialFile("myMaterials.material");

    vis->addResourceDirectory(vis->getResourceDir() + "/material/skybox/violentdays");
    vis->loadMaterialFile("violentdays.material");

    /// shdow setting
    vis->getSceneManager()->setShadowTechnique(Ogre::SHADOWTYPE_TEXTURE_ADDITIVE);
    vis->getSceneManager()->setShadowTextureSettings(2048, 3);

    /// scale related settings!! Please adapt it depending on your map size
    // beyond this distance, shadow disappears
    vis->getSceneManager()->setShadowFarDistance(10);
    // size of contact points and contact forces
    vis->setContactVisObjectSize(0.03, 0.6);
    // speed of camera motion in freelook mode
    vis->getCameraMan()->setTopSpeed(5);
}

void disturbance(std::vector<raisim::ArticulatedSystem *> A1,std::map<std::string, raisim::VisualObject> &objList,size_t start,size_t stop, size_t ctrlTick){
    Eigen::VectorXd pos = Eigen::MatrixXd::Zero(TOTAL_DOF+1,1);
    Eigen::VectorXd vel = Eigen::MatrixXd::Zero(TOTAL_DOF,1);
    A1.back()->getState(pos, vel);
    
    // ====================== External Forces ====================== //
    raisim::Vec<3> extForce;
    raisim::Mat<3,3> rot;
    raisim::Vec<3> dir;
    if (ctrlTick>=start && ctrlTick<stop){
        extForce = {0,-20,0}; // Pulse
        // extForce = {50*sin(4*ctrlTick*simfreq_raisim),0,0}; // Fwd Sine
        //extForce = {0,20*sin(4*ctrlTick*simfreq_raisim),0}; // Lat Sine
    } else{
        extForce = {0,0,0};
    }
    A1.back()->setExternalForce(1,extForce);
    dir = extForce;
    dir /= dir.norm();
    raisim::zaxisToRotMat(dir, rot);
    objList["extForceArrow"].offset = {pos(0),pos(1),pos(2)};
    objList["extForceArrow"].scale = {0.2,0.2,0.01*extForce.norm()};
    objList["extForceArrow"].rotationOffset = rot;
}

void plotGRFs(std::map<std::string, raisim::VisualObject>& objList, const std::vector<double>& GRFs, const std::vector<double>& pf, const std::vector<double>& contact, int agent_idx) {
    // Ensure the vectors are of the correct size
    if (GRFs.size() < 12 || pf.size() < 12 || contact.size() < 4) {
        std::cerr << "Error: Input vectors are of incorrect size." << std::endl;
        return;
    }

    for (int i = 0; i < 4; ++i) {
        raisim::Vec<3> dir;
        for (int j = 0; j < 3; ++j) {
            dir[j] = GRFs[3 * i + j];
        }

        // Normalize the direction vector if it is not a zero vector
        double norm = dir.norm();
        if (norm > 1e-6) {  // Check if the vector is non-zero to avoid division by zero
            dir /= norm;
        } else {
            dir.setZero();
        }

        // Convert direction vector to a rotation matrix that aligns the z-axis with the direction vector
        raisim::Mat<3, 3> rot;
        if (contact[i] == 1 && norm > 1e-6) {
            raisim::zaxisToRotMat(dir, rot);
        } else {
            rot.setIdentity();  // Set rotation to identity if no contact or zero norm
        }

        // Visual object key
        std::string objKey = "GRF" + std::to_string(agent_idx+1) + std::to_string(i + 1);
        objList[objKey].offset = {pf[3 * i], pf[3 * i + 1], pf[3 * i + 2]};
        objList[objKey].scale = {0.2, 0.2, 0.005 * norm};  // Scaling based on the norm of the GRF vector

        // Set the rotation matrix
        objList[objKey].rotationOffset = rot;
    }
}

void controller(size_t agent_idx, std::vector<raisim::ArticulatedSystem *> A1, std::vector<LocoWrapper*> loco_objs, std::vector<NMPC*> nmpc_objs, size_t controlTick, std::vector<int>& runMPC) {
    /////////////////////////////////////////////////////////////////////
    //////////////////////////// INITIALIZE
    /////////////////////////////////////////////////////////////////////
    size_t loco_kind = TROT;                        // Gait pattern to use
    size_t pose_kind = POSE_COMB;                   // Pose type to use (if loco_kind is set to POSE)
    size_t settling = 0.2 * ctrlHz;                 // Settling down
    size_t duration = 0.8 * ctrlHz;                 // Stand up
    size_t loco_start = settling + duration;        // Start the locomotion pattern

    double* tau;
    Eigen::VectorXd jointTorqueFF = Eigen::MatrixXd::Zero(TOTAL_DOF, 1);
    Eigen::VectorXd jointPosTotal = Eigen::MatrixXd::Zero(TOTAL_DOF + 1, 1); // +1 is for 4th Component of Quaternion 
    Eigen::VectorXd jointVelTotal = Eigen::MatrixXd::Zero(TOTAL_DOF, 1);
    static int conIndDes[4] = { 1 };

    /////////////////////////////////////////////////////////////////////
    //////////////////////////// UPDATE STATE
    /////////////////////////////////////////////////////////////////////

    A1[agent_idx]->getState(jointPosTotal, jointVelTotal);
    raisim::Mat<3, 3> rotMat;
    Eigen::Matrix3d rotE;
    A1[agent_idx]->getBaseOrientation(rotMat);
    double rotMatrixDouble[9];
    for (size_t j = 0; j < 9; ++j) {
        rotMatrixDouble[j] = rotMat[j];
        rotE(j) = rotMat[j];
    }
    jointVelTotal.segment(3, 3) = rotE.transpose() * jointVelTotal.segment(3, 3); // convert to body frame, like robot measurements

    double jpos[18], jvel[18];
    Eigen::Matrix<double, 3, 1> eul;
    Eigen::Matrix<double, 4, 1> quat;
    quat = jointPosTotal.block(3, 0, 4, 1);
    quat_to_XYZ(quat, eul);
    for (size_t k = 0; k < 3; ++k) {
        jpos[k] = jointPosTotal(k);
        jvel[k] = jointVelTotal(k);
        jpos[k + 3] = eul(k);
        jvel[k + 3] = jointVelTotal(k + 3);
    }
    for (size_t k = 6; k < 18; ++k) {
        jpos[k] = jointPosTotal(k + 1);
        jvel[k] = jointVelTotal(k);
    }

    int force[4] = { 0 };
    for (auto& con : A1[agent_idx]->getContacts()) {
        int conInd = con.getlocalBodyIndex();
        force[conInd / 3 - 1] = 0;
    }

    // if (agent_idx == 0) {
    //     kinEst0(force,conIndDes,jpos,jvel,rotE);
    // }
    // else
    // {
    //     kinEst1(force,conIndDes,jpos,jvel,rotE);
    // }
    /////////////////////////////////////////////////////////////////////
    //////////////////////////// CONTROL
    /////////////////////////////////////////////////////////////////////
    // Update the desired torques
    if (controlTick < settling) { // Settle down
        double temp[18] = { 0 };
        tau = temp;
        loco_objs[agent_idx]->initStandVars(jointPosTotal.block(0, 0, 3, 1), (int)duration);
    }
    else if (controlTick >= settling && controlTick <= loco_start) { // Start standing

        //================================================================================================
        //========================================== HIGH LEVEL ==========================================
        //================================================================================================

        updateData(GET_DATA, HL_DATA, &agentData[agent_idx], agent_idx);
        nmpc_objs[agent_idx]->updateState(agentData[agent_idx].q, 
                                        agentData[agent_idx].dq, 
                                        agentData[agent_idx].ind, 
                                        agentData[agent_idx].toePos, 
                                        agentData[1 - agent_idx].last_state[1 - agent_idx], 
                                        agentData[1 - agent_idx].GRFs_sol_.row(1 - agent_idx), 
                                        agentData[agent_idx].domain, 
                                        controlTick);

        if (runMPC[agent_idx] == 1 && agentData[agent_idx].domain < TOTALSTEPNUM)
        {
            runMPC[agent_idx] = 0;
            nmpc_objs[agent_idx]->run_NMPC();

            agentData[agent_idx].MPC_data_available = 1;
            agentData[agent_idx].cbf_active = nmpc_objs[agent_idx]->getCBF_status();
            agentData[agent_idx].alpha_COM = nmpc_objs[agent_idx]->get_alphaCOM();
            agentData[agent_idx].alpha_GRF = nmpc_objs[agent_idx]->get_alphaGRF();
            agentData[agent_idx].MPC_sol_ = nmpc_objs[agent_idx]->get_MPCsol();
            agentData[agent_idx].GRFs_sol_.row(agent_idx) = nmpc_objs[agent_idx]->get_GRFsol();
            agentData[agent_idx].last_state[agent_idx] = nmpc_objs[agent_idx]->get_lastState();
            // std::cout << "GRFs_sol_ (A1_Sim): agent (" << agent_idx << "): " << agentData[agent_idx].GRFs_sol_.row(agent_idx) << std::endl;
            // std::cin.get();
            // std::cout << "RUN MPC while in STAND" << std::endl;
        }
        updateData(SET_DATA, HL_DATA, &agentData[agent_idx], agent_idx);

        //================================================================================================
        //========================================== LOW LEVEL ===========================================
        //================================================================================================

        updateData(GET_DATA, LL_DATA, &agentData[agent_idx], agent_idx);
        loco_objs[agent_idx]->set_MPC_DATA( agentData[agent_idx].alpha_COM,
                                            agentData[agent_idx].alpha_GRF, 
                                            agentData[agent_idx].MPC_sol_, 
                                            agentData[agent_idx].MPC_data_available,
                                            agentData[agent_idx].cbf_active);
        
        loco_objs[agent_idx]->calcTau(jpos, jvel, rotMatrixDouble, force, STAND, controlTick, &runMPC[agent_idx]);
        tau = loco_objs[agent_idx]->getTorque();

        agentData[agent_idx].MPC_data_available = loco_objs[agent_idx]->get_MPC_data_available();
        const int* contactMat = loco_objs[agent_idx]->getConDes();
        Eigen::Matrix<double, 3, 4> toePosTmp = loco_objs[agent_idx]->get_toePos();
        agentData[agent_idx].toePos = toePosTmp;
        agentData[agent_idx].domain = loco_objs[agent_idx]->getDomain();
        agentData[agent_idx].domain_change = loco_objs[agent_idx]->hasDomainChanged();
        // std::cout << "agentData" << agent_idx << ".domain: " << agentData[agent_idx].domain << std::endl;
        memcpy(agentData[agent_idx].ind, contactMat, 4 * sizeof(int));
        updateData(SET_DATA, LL_DATA, &agentData[agent_idx], agent_idx);
    }
    else if (controlTick > loco_start) { // Start locomotion

        //================================================================================================
        //========================================== HIGH LEVEL ==========================================
        //================================================================================================
        updateData(GET_DATA, HL_DATA, &agentData[agent_idx], agent_idx);
        nmpc_objs[agent_idx]->updateState(agentData[agent_idx].q, 
                                  agentData[agent_idx].dq, 
                                  agentData[agent_idx].ind, 
                                  agentData[agent_idx].toePos, 
                                  agentData[1 - agent_idx].last_state[1 - agent_idx], 
                                  agentData[1 - agent_idx].GRFs_sol_.row(1 - agent_idx), 
                                  agentData[agent_idx].domain, 
                                  controlTick);
        // std::cout << "Debugging GRFs_sol_: agent (" << 1 - agent_idx << "): " << nmpc_objs[agent_idx]->get_GRFsol() << std::endl;

        if ((runMPC[agent_idx] == 1 || sec_time[agent_idx]) && agentData[agent_idx].domain < TOTALSTEPNUM)
        {
            sec_time[agent_idx] = 0;
            runMPC[agent_idx] = 0;

            nmpc_objs[agent_idx]->run_NMPC();
            agentData[agent_idx].MPC_data_available = 1;
            agentData[agent_idx].cbf_active = nmpc_objs[agent_idx]->getCBF_status();
            agentData[agent_idx].alpha_COM = nmpc_objs[agent_idx]->get_alphaCOM();
            agentData[agent_idx].alpha_GRF = nmpc_objs[agent_idx]->get_alphaGRF();

            // std::cout << "alpha_GRF_:\n" << agentData[agent_idx].alpha_GRF << std::endl; 

            agentData[agent_idx].MPC_sol_ = nmpc_objs[agent_idx]->get_MPCsol();
            agentData[agent_idx].GRFs_sol_.row(agent_idx) = nmpc_objs[agent_idx]->get_GRFsol();
            agentData[agent_idx].last_state[agent_idx] = nmpc_objs[agent_idx]->get_lastState();
        }

        updateData(SET_DATA, HL_DATA, &agentData[agent_idx], agent_idx);

        //================================================================================================
        //========================================== LOW LEVEL ==========================================
        //================================================================================================

        updateData(GET_DATA, LL_DATA, &agentData[agent_idx], agent_idx);
        loco_objs[agent_idx]->set_MPC_DATA( agentData[agent_idx].alpha_COM,
                                            agentData[agent_idx].alpha_GRF, 
                                            agentData[agent_idx].MPC_sol_, 
                                            agentData[agent_idx].MPC_data_available, 
                                            agentData[agent_idx].cbf_active);

        loco_objs[agent_idx]->calcTau(jpos, jvel, rotMatrixDouble, force, loco_kind, controlTick, &runMPC[agent_idx]);
        tau = loco_objs[agent_idx]->getTorque();

        agentData[agent_idx].MPC_data_available = loco_objs[agent_idx]->get_MPC_data_available();
        const int* contactMat = loco_objs[agent_idx]->getConDes();
        Eigen::Matrix<double, 3, 4> toePosTmp = loco_objs[agent_idx]->get_toePos();
        agentData[agent_idx].toePos = toePosTmp;
        agentData[agent_idx].domain = loco_objs[agent_idx]->getDomain();
        agentData[agent_idx].domain_change = loco_objs[agent_idx]->hasDomainChanged();
        memcpy(agentData[agent_idx].ind, contactMat, 4 * sizeof(int));
        updateData(SET_DATA, LL_DATA, &agentData[agent_idx], agent_idx);
    }

    memcpy(agentData[agent_idx].q, jpos, 18 * sizeof(double));
    memcpy(agentData[agent_idx].dq, jvel, 18 * sizeof(double));

    const int* contactMat = loco_objs[agent_idx]->getConDes();
    memcpy(agentData[agent_idx].ind, contactMat, 4 * sizeof(int));
    // memcpy(conIndDes, contactMat, 4 * sizeof(int));

    jointTorqueFF = Eigen::Map<Eigen::Matrix<double, 18, 1>>(tau, 18);
    jointTorqueFF.block(0, 0, 6, 1).setZero();

    // Set the desired torques
    A1[agent_idx]->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
    A1[agent_idx]->setGeneralizedForce(jointTorqueFF);
}

void init_vars()
{
    HLData.resize(NUMBER_OF_AGENTS);
    LLData.resize(NUMBER_OF_AGENTS);
    HLData_Backup.resize(NUMBER_OF_AGENTS);
    LLData_Backup.resize(NUMBER_OF_AGENTS);
}

Eigen::MatrixXd loadObstacles(std::string key) {
    YAML::Node config = YAML::LoadFile("../src/config.yaml");
    std::vector<std::pair<double, double>> obstaclesVec;

    for (const auto& node : config[key]) {
        double x = node[0].as<double>();
        double y = node[1].as<double>();
        obstaclesVec.emplace_back(x, y);
    }

    // Convert vector to Eigen::MatrixXd
    Eigen::MatrixXd obstacles(2, obstaclesVec.size());
    for (size_t i = 0; i < obstaclesVec.size(); ++i) {
        obstacles(0, i) = obstaclesVec[i].first;
        obstacles(1, i) = obstaclesVec[i].second;
    }

    return obstacles;

}

int main(int argc, char *argv[]) {
    init_vars();
    const int NUMBER_OF_SIMS = 1;
    const float threshold = 0.4;

    // ============================================================ //
    // =================== SETUP RAISIM/VISUALS =================== //
    // ============================================================ //

    raisim::World::setActivationKey(raisim::loadResource("activation.raisim"));
    raisim::World world;
    world.setTimeStep(simfreq_raisim);

    raisim::OgreVis *vis = raisim::OgreVis::get();

    /// these method must be called before initApp
    vis->setWorld(&world);
    vis->setWindowSize(1792, 1200); // Should be evenly divisible by 16!!
    vis->setImguiSetupCallback(imguiSetupCallback); // These 2 lines make the interactable gui visible
    vis->setImguiRenderCallback(imguiRenderCallBack);
    vis->setKeyboardCallback(raisimKeyboardCallback);
    vis->setSetUpCallback(setupCallback);
    vis->setAntiAliasing(2);

    /// starts visualizer thread
    vis->initApp();
    /// create raisim objects
    raisim::TerrainProperties terrainProperties;
    terrainProperties.frequency = 0.0;
    terrainProperties.zScale = 0.0;
    terrainProperties.xSize = 300.0;
    terrainProperties.ySize = 300.0;
    terrainProperties.xSamples = 50;
    terrainProperties.ySamples = 50;
    terrainProperties.fractalOctaves = 0;
    terrainProperties.fractalLacunarity = 0.0;
    terrainProperties.fractalGain = 0.0;

    raisim::HeightMap *ground = world.addHeightMap(0.0, 0.0, terrainProperties);
    vis->createGraphicalObject(ground, "terrain", "checkerboard_blue");
    world.setDefaultMaterial(0.8, 0.0, 0.0); //surface friction could be 0.8 or 1.0
    vis->addVisualObject("extForceArrow", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);

    // Add path visualization for both agents 
    // Define the reference path points (example points)
    
    // Create visual object for the reference path

    // Define the reference path points (example points)
    // std::vector<Eigen::Vector3d> referencePath = {
    //     {0.0, 0.0, 0.0},
    //     {1.0, 0.5, 0.0},
    //     {2.0, 1.0, 0.0},
    //     {3.0, 0.5, 0.0},
    //     {4.0, 0.0, 0.0}
    // };

    // // Create and configure cylinders for each segment of the path
    // for (size_t i = 1; i < referencePath.size(); ++i) {
    //     Eigen::Vector3d start = referencePath[i - 1];
    //     Eigen::Vector3d end = referencePath[i];
    //     Eigen::Vector3d midPoint = (start + end) / 2.0;
    //     Eigen::Vector3d direction = (end - start).normalized();
    //     double length = (end - start).norm();

    //     // Create orientation quaternion for the cylinder
    //     Eigen::Quaterniond orientation = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), direction);

    //     // Create a cylinder object in the world
    //     auto* cylinder = world.addCylinder(0.05, length, 0.1);

    //     // Set position and orientation of the cylinder
    //     cylinder->setPosition(midPoint.x(), midPoint.y(), midPoint.z());
    //     cylinder->setOrientation(orientation.w(), orientation.x(), orientation.y(), orientation.z());

    //     // Create a visual object for the cylinder
    //     vis->createGraphicalObject(cylinder, "pathref" + std::to_string(i), "yellow");
    // }



    // Foot force visualization arrows
    vis->addVisualObject("GRF11", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF12", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF13", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF14", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF21", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF22", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF23", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF24", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);

    auto& list = vis->getVisualObjectList();

    float box_width = 0.15;
    float box_z = 0.15;
    auto boxd = world.addBox(box_width,box_width,0.3,0.1);
    vis->createGraphicalObject(boxd, "boxd", "red");

    // Mimic the rough terrain

    std::vector<std::string> colors = {"red", "green", "blue", "yellow", "orange", "purple", "pink", "brown", "grey", "white"};
    std::random_device rd1, rd2, rd3;  
    std::mt19937 gen1(rd1()), gen2(rd2()), gen3(rd3());
    std::uniform_real_distribution<> dis_height(0.0, 0.02); // Random heights between 0 and 1 cm
    std::uniform_real_distribution<> dis_width(0.03, 0.2); // Random width between 3 and 7 cm
    std::uniform_real_distribution<> dis_length(1, 5); // Random length between 3 and 7 cm
    std::uniform_int_distribution<> dis_color(0, colors.size() - 1); // Distribution for color indices

    int num_blocks = 24; // Number of blocks
    std::vector<raisim::Box*> wood_blocks;

    int num_rows = num_blocks/1; // Number of rows for block placement
    float spacing_factor = 1.5; // Adjust spacing between blocks

    if(rough_terrain_en) {
        for (size_t i = 0; i < num_blocks; i++) {
            float wood_height = dis_height(gen1); // Generate a random height for each box
            float wood_width = dis_width(gen2);   // Generate a random width for each box
            float wood_length = dis_length(gen3); // Generate a random length for each box
            int color_index = dis_color(gen1);    // Generate a random color index

            wood_blocks.push_back(world.addBox(wood_width, wood_length, wood_height, 0.1)); // Add the box to the world with random dimensions
            vis->createGraphicalObject(wood_blocks.back(), "wood_block" + std::to_string(i), colors[color_index]);

            // Set the position of the box in a grid pattern
            int row = i / num_rows;
            int col = i % num_rows;
            wood_blocks.back()->setPosition(col * wood_width * spacing_factor, row * wood_length * spacing_factor, wood_height / 2);
        }
    };

    // ============================================================ //
    // ===================== SETUP Obstacles ====================== //
    // ============================================================ //

    std::vector<raisim::Cylinder*> obstacles;
    if (obstacle_en)
    {
        for (size_t i = 0; i < NUMBER_OF_OBS; i++)
        {
            obstacles.push_back(world.addCylinder(box_width,box_width,0.3));
            vis->createGraphicalObject(obstacles.back(), "box"+std::to_string(i), "red");

            // Visualize non-real obstacles
            std::string objName = "cylinder" + std::to_string(i + 1); // Unique name for each cylinder
            vis->addVisualObject(objName, "cylinderMesh", "yellow", {0.15, 0.15, 0.2}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);  
        }
    }

    // ============================================================ //
    // ======================= SETUP Robot ======================== //
    // ============================================================ //

    Eigen::Matrix<double, 2*NUMBER_OF_AGENTS, 1> Pstart;
    Pstart << 0.0, 0.0, 0.0, -1.0;

    std::vector<raisim::ArticulatedSystem*> A1;

    for (size_t agent_idx = 0; agent_idx < NUMBER_OF_AGENTS; ++agent_idx) {

        A1.push_back(world.addArticulatedSystem(raisim::loadResource("A1/A1_modified_new.urdf")));
        vis->createGraphicalObject(A1.back(), "A1_" + std::to_string(agent_idx));
        A1.back()->setName("A1_Robot_" + std::to_string(agent_idx));

        // Set the generalized coordinates for the current agent
        double x_pos = Pstart(2 * agent_idx);
        double y_pos = Pstart(2 * agent_idx + 1);
        A1.back()->setGeneralizedCoordinate({x_pos, y_pos, 0.12, 1, 0, 0, 0, 0.0, M_PI / 3, -2.6, 0.0, M_PI / 3, -2.6, 0.0, M_PI / 3, -2.6, 0.0, M_PI / 3, -2.6});
        // A1.back()->setGeneralizedCoordinate({Pstart(0), Pstart(1), 0.12, 1, 0, 0, 0,0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6});
        A1.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);
    }
    
    // ============================================================ //
    // ========================= SIM EPISODE ====================== //
    // ============================================================ //

    for (size_t sim = 0; sim < NUMBER_OF_SIMS; ++sim)
    {
        boxd->setPosition(250, 0.1,box_z);
        // Define the matrices
        Eigen::Matrix<double, 2, NUMBER_OF_OBS> Pobs;
        Eigen::Matrix<double, 2, NUMBER_OF_OBS> Pobs_real;


        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        
        // Uniform distributions for the required ranges
        std::uniform_real_distribution<> dis_x(1, 9.0); // For x values 
        std::uniform_real_distribution<> dis_y(-3.0, 3.0); // For y values 
        std::uniform_real_distribution<> dis_uncertainty(-0.2, 0.2); // For uncertainties

        
        srand(static_cast<unsigned int>(time(nullptr)));

        for (int col = 0; col < Pobs.cols(); ++col) {
            // Generate random values in the specified ranges
            double randX = dis_x(gen);
            double randY = dis_y(gen);

            // Assign these values to Pobs
            Pobs(0, col) = randX;
            Pobs(1, col) = randY;

            if(!obstacle_en)
            {
                // Add 100 in y values to move obstacles out of the way 
                Pobs.row(1).array() += 100;
            }

            // Generate uncertainty
            double uncertaintyX = dis_uncertainty(gen);
            double uncertaintyY = dis_uncertainty(gen);

            // Add uncertainty to create Pobs_real
            Pobs_real(0, col) = Pobs(0, col) + uncertaintyX;
            Pobs_real(1, col) = Pobs(1, col) + uncertaintyY;

            // Ensure Pobs_real stays within bounds
            Pobs_real(0, col) = std::min(std::max(Pobs_real(0, col), 1.0), 9.0);
            Pobs_real(1, col) = std::min(std::max(Pobs_real(1, col), -3.0), 3.0);
        }
        // Pobs <<  2.5,   2.5,     4.5,   4.5,     6.5,   6.5,     8.5,   8.5,
        //             1,      0,       1,      0,       1,      0,       1,      0;
        
        std::cout << "Loading Pobs from YAML...\n";
        Pobs = loadObstacles("Pobs");
        std::cout << "Loading Pobs_real from YAML...\n";
        Pobs_real = loadObstacles("Pobs_real");


    
        // ============================================================ //
        // =================== SETUP ENVIRONMENT ====================== //
        // ============================================================ //
        
        std::cout << "Pobs:\n" << Pobs << "\n\n";
        std::cout << "Pobs_real (with uncertainty):\n" << Pobs_real << "\n";

    
        // Populate the obstacles 
        int obs_iter = 0;

        auto& list_nonreal_obs = vis->getVisualObjectList();

        if (obstacle_en)
        {
            for (auto it = obstacles.begin(); it != obstacles.end(); ++it) {
                raisim::Cylinder* obstacle = *it;
                // Use 'obstacle' as needed
                obstacle->setPosition(Pobs_real(0,obs_iter),Pobs_real(1,obs_iter),box_z);
                // obstacle->setPosition(Pobs(0,obs_iter),Pobs(1,obs_iter),box_z);
                obstacle->setMass(0.2);
                obstacle->setOrientation(0,0,0,1);
                obs_iter++;
            }
            
            for (int i = 0; i < Pobs.cols(); ++i) {
                std::string objKey = "cylinder" + std::to_string(i + 1);
                if (list_nonreal_obs.find(objKey) != list_nonreal_obs.end()) {
                    double x = Pobs(0, i); // x coordinate of the i-th obstacle
                    double y = Pobs(1, i); // y coordinate of the i-th obstacle
                    list_nonreal_obs[objKey].offset = {x, y, 0.0}; // Set the position
                }
            }
        }

        std::vector<LocoWrapper*> loco_objs(NUMBER_OF_AGENTS);
        std::vector<NMPC*> nmpc_objs(NUMBER_OF_AGENTS);
        
        for (size_t agent_idx = 0; agent_idx < NUMBER_OF_AGENTS; ++agent_idx) {

            // Initialize the LocoWrapper for each agent
            loco_objs[agent_idx] = new LocoWrapper(argc, argv);
            loco_objs[agent_idx]->setAgentID(agent_idx); // Set the agent ID
            loco_objs[agent_idx]->setPstart(Pstart);
            loco_objs[agent_idx]->setPobs(Pobs);

            // Initialize the NMPC for each agent
            nmpc_objs[agent_idx] = new NMPC();
            nmpc_objs[agent_idx]->init();
            nmpc_objs[agent_idx]->setAgentID(agent_idx); // Set the agent ID
            nmpc_objs[agent_idx]->setPstart(Pstart);
            nmpc_objs[agent_idx]->setPobs(Pobs);
            nmpc_objs[agent_idx]->setPobs_real(Pobs_real);
            nmpc_objs[agent_idx]->generateReferenceTrajectory();
        }

        // ============================================================ //
        // ================= VIEW AND RECORDING OPTIONS =============== //
        // ============================================================ //
        raisim::gui::showContacts = false;
        raisim::gui::showForces = false;
        raisim::gui::showCollision = false;
        raisim::gui::showBodies = true;

        std::string cameraview = "top";
        bool panX = true;                // Pan view with robot during walking (X direction)
        bool panY = true;                // Pan view with robot during walking (Y direction)
        bool record = true;             // Record?
        double startTime = 0*ctrlHz;    // RecoPstartrding start time
        // double simlength = nmpc->params.sim_length/nmpc->params.dt; //60*ctrlHz;   // Sim end time
        double simlength = SIM_TIME*ctrlHz;
        double fps = 30;            
        std::string directory = "~/Desktop/Raisim_Simulations/";
        // std::string filename = "Payload_Inplace";
        std::string filename = "SRB_NMPC";
        // std::string filename = "inplace_sim";

        // ============================================================ //
        // ========================= VIEW SETUP ======================= //
        // ============================================================ //
        if(cameraview == "iso"){
            vis->getCameraMan()->getCamera()->setPosition(-1, -3, 0.5);
            vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(4.5*Pi/6-Pi/2));
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
        }else if(cameraview == "isoside"){
            vis->getCameraMan()->getCamera()->setPosition(1.1, -2, 0.5);
            vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(4*Pi/6-Pi/2));
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
        }else if(cameraview == "side"){
            vis->getCameraMan()->getCamera()->setPosition(0, -2, 0.5);
            vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(0));
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
        }else if(cameraview == "front"){
            vis->getCameraMan()->getCamera()->setPosition(2, 0, 0.5);
            vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(Pi/2));
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
        }else if(cameraview == "top"){
            vis->getCameraMan()->getCamera()->setPosition(2, -1, 6);
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(0));
        }else{
            vis->getCameraMan()->getCamera()->setPosition(1, -3, 2.5);
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(1.0));
        }
        unsigned long mask = 0;
        if(raisim::gui::showBodies) mask |= raisim::OgreVis::RAISIM_OBJECT_GROUP;
        if(raisim::gui::showCollision) mask |= raisim::OgreVis::RAISIM_COLLISION_BODY_GROUP;
        if(raisim::gui::showContacts) mask |= raisim::OgreVis::RAISIM_CONTACT_POINT_GROUP;
        if(raisim::gui::showForces) mask |= raisim::OgreVis::RAISIM_CONTACT_FORCE_GROUP;
        vis->setVisibilityMask(mask);

        if(panX) raisim::gui::panViewX = panX;
        if(panY) raisim::gui::panViewY = panY;
        
        // ============================================================ //
        // ========================== RUN SIM ========================= //
        // ============================================================ //
        const std::string name = directory+filename+"_"+cameraview+".mp4";
        vis->setDesiredFPS(fps);
        long simcounter = 0;
        static bool added = false;

        while (!vis->getRoot()->endRenderingQueued() && simcounter < simlength-1){

            for (int agent_idx = 0; agent_idx < NUMBER_OF_AGENTS; agent_idx++)
            {
                controller(agent_idx, A1, loco_objs, nmpc_objs, simcounter, runMPC);
                std::vector<double> GRF = nmpc_objs[agent_idx]->get_GRF();
                std::vector<double> feet_vec = nmpc_objs[agent_idx]->get_Pf();
                std::vector<double> contacts = nmpc_objs[agent_idx]->get_contacts();
                plotGRFs(list, GRF, feet_vec, contacts, agent_idx);
            }

            // MPC_Agent0->updateDistance_to_fail();
            // MPC_Agent1->updateDistance_to_fail();

            world.integrate();        
            
            if (simcounter%30 == 0)
                vis->renderOneFrame();
            
            if (!vis->isRecording() & record & simcounter>=startTime)
                vis->startRecordingVideo(name);
            
            auto currentPos = vis->getCameraMan()->getCamera()->getPosition();
            if (raisim::gui::panViewX){
                Eigen::VectorXd jointPosTotal(18 + 1);
                Eigen::VectorXd jointVelTotal(18);
                jointPosTotal.setZero();
                jointVelTotal.setZero();
                A1.back()->getState(jointPosTotal, jointVelTotal);
                if (cameraview=="front"){
                    currentPos[0] = jointPosTotal(0)+2;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="side"){
                    currentPos[0] = jointPosTotal(0);
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="iso"){
                    currentPos[0] = jointPosTotal(0)+2;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="isoside"){
                    currentPos[0] = jointPosTotal(0)+1.1;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                }
            }
            if (raisim::gui::panViewY){
                Eigen::VectorXd jointPosTotal(18 + 1);
                Eigen::VectorXd jointVelTotal(18);
                jointPosTotal.setZero();
                jointVelTotal.setZero();
                A1.back()->getState(jointPosTotal, jointVelTotal);
                if (cameraview=="front"){
                    currentPos[1] = jointPosTotal(1);
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="side"){
                    currentPos[1] = jointPosTotal(1)-2;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="iso"){
                    currentPos[1] = jointPosTotal(1)-1;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="isoside"){
                    currentPos[1] = jointPosTotal(1)-2;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                }
            }
            simcounter++;
        }
        
        // for (size_t idx = 0; idx < NUMBER_OF_AGENTS; idx++)
        // {
        //     nmpc_objs[idx]->logData();
        // }

        // outFile << failingDistance0 << "," << failingDistance1 << std::endl;

        // Cleanup
        for (size_t agent_idx = 0; agent_idx < NUMBER_OF_AGENTS; ++agent_idx) {
            delete A1[agent_idx];
            delete loco_objs[agent_idx];
            delete nmpc_objs[agent_idx];
        }
        
    }
    // outFile.close();

    // End recording if still recording
    if (vis->isRecording())
        vis->stopRecordingVideoAndSave();
    
    /// terminate the app
    vis->closeApp();
    return 0;
}