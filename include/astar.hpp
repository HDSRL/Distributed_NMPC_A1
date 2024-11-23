#include <queue>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

struct Node {
    int x, y;
    double g, h;
    Node* parent;

    Node(int x, int y, double g, double h, Node* parent)
        : x(x), y(y), g(g), h(h), parent(parent) {}

    double f() const { return g + h; }
};

struct CompareNode {
    bool operator()(const Node* a, const Node* b) const {
        return a->f() > b->f();
    }
};

double heuristic(int x1, int y1, int x2, int y2) {
    return std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));  // Euclidean distance
}

bool isValid(int x, int y, const Eigen::MatrixXd& obstacle_map) {
    if (x < 0 || x >= obstacle_map.rows() || y < 0 || y >= obstacle_map.cols()) {
        return false;
    }
    return obstacle_map(x, y) == 0;
}

std::vector<Node*> get_neighbors(Node* node, const Eigen::MatrixXd& obstacle_map, int goal_x, int goal_y) {
    std::vector<Node*> neighbors;
    std::vector<std::pair<int, int>> directions = {
        {1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {-1, -1}, {1, -1}, {-1, 1}  // Include diagonal directions
    };
    for (auto dir : directions) {
        int nx = node->x + dir.first;
        int ny = node->y + dir.second;
        if (isValid(nx, ny, obstacle_map)) {
            double g = node->g + std::hypot(dir.first, dir.second); // Use Euclidean distance for g-cost
            double h = heuristic(nx, ny, goal_x, goal_y);
            neighbors.push_back(new Node(nx, ny, g, h, node));
        }
    }
    return neighbors;
}

std::vector<Eigen::Vector2d> a_star(int start_x, int start_y, int goal_x, int goal_y, const Eigen::MatrixXd& obstacle_map) {
    std::priority_queue<Node*, std::vector<Node*>, CompareNode> open_list;
    std::unordered_map<int, Node*> all_nodes;
    Node* start_node = new Node(start_x, start_y, 0, heuristic(start_x, start_y, goal_x, goal_y), nullptr);
    open_list.push(start_node);
    all_nodes[start_x * obstacle_map.cols() + start_y] = start_node;

    while (!open_list.empty()) {
        Node* current = open_list.top();
        open_list.pop();

        if (current->x == goal_x && current->y == goal_y) {
            std::vector<Eigen::Vector2d> path;
            while (current != nullptr) {
                path.push_back(Eigen::Vector2d(current->x, current->y));
                current = current->parent;
            }
            std::reverse(path.begin(), path.end());
            return path;
        }

        for (Node* neighbor : get_neighbors(current, obstacle_map, goal_x, goal_y)) {
            int hash = neighbor->x * obstacle_map.cols() + neighbor->y;
            if (all_nodes.find(hash) == all_nodes.end() || neighbor->g < all_nodes[hash]->g) {
                open_list.push(neighbor);
                all_nodes[hash] = neighbor;
            }
        }
    }

    return std::vector<Eigen::Vector2d>(); // Return empty path if no path found
}

std::vector<Eigen::Vector2d> filter_waypoints(const std::vector<Eigen::Vector2d>& path) {
    if (path.size() <= 2) return path; // No filtering needed for paths with 2 or fewer points

    std::vector<Eigen::Vector2d> filtered_path;
    filtered_path.push_back(path.front());

    Eigen::Vector2d prev_direction = path[1] - path[0];
    for (size_t i = 2; i < path.size(); ++i) {
        Eigen::Vector2d current_direction = path[i] - path[i - 1];
        if (current_direction.normalized() != prev_direction.normalized()) {
            filtered_path.push_back(path[i - 1]);
        }
        prev_direction = current_direction;
    }

    filtered_path.push_back(path.back());
    return filtered_path;
}
