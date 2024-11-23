#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <random>

struct RRTNode {
    Eigen::Vector2d position;
    RRTNode* parent;
};

class RRT {
public:
    RRT(const Eigen::Vector2d& start, const Eigen::Vector2d& goal, const Eigen::MatrixXd& obstacles, double step_size, double goal_bias, int max_iter)
        : start(start), goal(goal), obstacles(obstacles), step_size(step_size), goal_bias(goal_bias), max_iter(max_iter) {}

    std::vector<Eigen::Vector2d> plan() {
        std::vector<RRTNode*> tree;
        RRTNode* root = new RRTNode{start, nullptr};
        tree.push_back(root);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);
        std::uniform_real_distribution<> x_dis(0, 100);
        std::uniform_real_distribution<> y_dis(0, 100);

        for (int i = 0; i < max_iter; ++i) {
            Eigen::Vector2d rand_point;
            if (dis(gen) < goal_bias) {
                rand_point = goal;
            } else {
                rand_point = Eigen::Vector2d(x_dis(gen), y_dis(gen));
            }

            RRTNode* nearest_node = nullptr;
            double nearest_dist = std::numeric_limits<double>::max();
            for (RRTNode* node : tree) {
                double dist = (node->position - rand_point).norm();
                if (dist < nearest_dist) {
                    nearest_node = node;
                    nearest_dist = dist;
                }
            }

            Eigen::Vector2d direction = (rand_point - nearest_node->position).normalized();
            Eigen::Vector2d new_point = nearest_node->position + step_size * direction;

            if (isCollisionFree(nearest_node->position, new_point)) {
                RRTNode* new_node = new RRTNode{new_point, nearest_node};
                tree.push_back(new_node);

                if ((new_point - goal).norm() < step_size) {
                    return extractPath(new_node);
                }
            }
        }

        return std::vector<Eigen::Vector2d>(); // Return empty path if not found
    }

private:
    bool isCollisionFree(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
        for (int i = 0; i < obstacles.cols(); ++i) {
            Eigen::Vector2d obs = obstacles.col(i);
            double dist = pointLineDistance(p1, p2, obs);
            if (dist < 1.0) {
                return false;
            }
        }
        return true;
    }

    double pointLineDistance(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p) {
        Eigen::Vector2d line = p2 - p1;
        Eigen::Vector2d point_to_line_start = p - p1;
        double t = point_to_line_start.dot(line) / line.squaredNorm();
        t = std::max(0.0, std::min(1.0, t)); // Clamp t to [0, 1]
        Eigen::Vector2d projection = p1 + t * line;
        return (p - projection).norm();
    }

    std::vector<Eigen::Vector2d> extractPath(RRTNode* node) {
        std::vector<Eigen::Vector2d> path;
        while (node != nullptr) {
            path.push_back(node->position);
            node = node->parent;
        }
        std::reverse(path.begin(), path.end());
        return path;
    }

    Eigen::Vector2d start;
    Eigen::Vector2d goal;
    Eigen::MatrixXd obstacles;
    double step_size;
    double goal_bias;
    int max_iter;
};
