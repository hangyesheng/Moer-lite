#include "BVH.h"

struct BVH::BVHNode {
    //* todo BVH节点结构设计
    BVHNode* left;
    BVHNode* right;
    AABB box;
    int firstShapeOffset;
    int nShape = 0;
    int splitAxis;
};
void BVH::build() {
    AABB sceneBox;
    for (const auto& shape : this->shapes) {
        //* 自行实现的加速结构请务必对每个shape调用该方法，以保证TriangleMesh构建内部加速结构
        //* 由于使用embree时，TriangleMesh::getAABB不会被调用，因此出于性能考虑我们不在TriangleMesh
        //* 的构造阶段计算其AABB，因此当我们将TriangleMesh的AABB计算放在TriangleMesh::initInternalAcceleration中
        //* 所以请确保在调用TriangleMesh::getAABB之前先调用TriangleMesh::initInternalAcceleration
        shape->initInternalAcceleration();
        sceneBox.Expand(shape->getAABB());
    }
    this->boundingBox = sceneBox;
    //* todo 完成BVH构建
    this->root = this->recursiveBuild(this->shapes, 0, this->shapes.size());
}

BVH::BVHNode* BVH::recursiveBuild(std::vector < std::shared_ptr < Shape>>& shapes, int l, int r) {
    if (l >= r) {
        return nullptr;
    }
    else if (r - l <= this->bvhLeafMaxSize) {
        AABB sceneBox;
        for (int i = l; i < r; i++) {
            sceneBox.Expand(shapes[i]->getAABB());
        }
        return new BVHNode{ nullptr, nullptr, sceneBox, l, r - l, 0 };
    }
    else {
        AABB sceneBox;
        std::vector<float> sum(3, 0);
        std::vector<float> sum_squared(3, 0);
        for (int i = l; i < r; i++) {
            AABB shapeBox = shapes[i]->getAABB();
            Point3f center = shapeBox.Center();
            for (int dimension = 0; dimension < 3; dimension++) {
                sum[dimension] += center[dimension];
                sum_squared[dimension] += center[dimension] * center[dimension];
            }
            sceneBox.Expand(shapeBox);
        }
        std::vector<float> variance(3, 0);
        for (int dimension = 0; dimension < 3; dimension++) {
            float mean = sum[dimension] / (r - l);
            variance[dimension] = sum_squared[dimension] / (r - l) - mean * mean;
        }
        int splitAxis = 0;
        if (variance[1] > variance[0] && variance[1] > variance[2]) {
            splitAxis = 1;
        }
        else if (variance[2] > variance[0] && variance[2] > variance[1]) {
            splitAxis = 2;
        }
        std::sort(shapes.begin() + l, shapes.begin() + r, [splitAxis](const std::shared_ptr<Shape>& a, const std::shared_ptr<Shape>& b) {
            return a->getAABB().Center()[splitAxis] < b->getAABB().Center()[splitAxis];
            });
        return new BVHNode{ recursiveBuild(shapes, l, (l + r) / 2), recursiveBuild(shapes, (l + r) / 2, r), sceneBox, l, r - l, splitAxis };
    }
}

bool BVH::rayIntersect(Ray& ray, int* geomID, int* primID, float* u, float* v) const {
    //* todo 完成BVH求交

    *geomID = -1;
    if (this->root != nullptr) {
        std::vector<BVHNode*> stack;
        stack.push_back(this->root);
        while (!stack.empty()) {
            BVHNode* node = stack.back();
            stack.pop_back();
            if (node->box.RayIntersect(ray)) {
                if (node->left == nullptr && node->right == nullptr) {
                    for (int i = node->firstShapeOffset; i < node->firstShapeOffset + node->nShape; i++) {
                        if (shapes[i]->rayIntersectShape(ray, primID, u, v)) {
                            *geomID = i;
                        }
                    }
                }
                else {
                    if (ray.direction[node->splitAxis] > 0) {
                        stack.push_back(node->right);
                        stack.push_back(node->left); // 先访问左节点
                    }
                    else {
                        stack.push_back(node->left);
                        stack.push_back(node->right); // 先访问右节点
                    }
                }
            }
        }
    }
    if (*geomID != -1) return true;

    return false;
}


