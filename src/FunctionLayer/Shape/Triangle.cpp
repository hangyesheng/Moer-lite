#include "Triangle.h"
#include <FunctionLayer/Acceleration/Linear.h>
//--- Triangle ---
Triangle::Triangle(int _primID, int _vtx0Idx, int _vtx1Idx, int _vtx2Idx,
  const TriangleMesh* _mesh)
  : primID(_primID), vtx0Idx(_vtx0Idx), vtx1Idx(_vtx1Idx), vtx2Idx(_vtx2Idx),
  mesh(_mesh) {
  Point3f vtx0 = mesh->transform.toWorld(mesh->meshData->vertexBuffer[vtx0Idx]),
    vtx1 = mesh->transform.toWorld(mesh->meshData->vertexBuffer[vtx1Idx]),
    vtx2 = mesh->transform.toWorld(mesh->meshData->vertexBuffer[vtx2Idx]);
  boundingBox.Expand(vtx0);
  boundingBox.Expand(vtx1);
  boundingBox.Expand(vtx2);
  this->geometryID = mesh->geometryID;
}

bool Triangle::rayIntersectShape(Ray& ray, int* primID, float* u,
  float* v) const {
  //* todo 实现三角形与光线求交

  // O+tD = (1-u-v)V0 + uV1 + vV2
  // [-D, V1-V0, V2-V0] [t,u,v]^T = O-V0 

  Point3f origin = ray.origin;
  Vector3f direction = ray.direction;

  Point3f v_0 = this->mesh->meshData->vertexBuffer[this->vtx0Idx];
  Point3f v_1 = this->mesh->meshData->vertexBuffer[this->vtx1Idx];
  Point3f v_2 = this->mesh->meshData->vertexBuffer[this->vtx2Idx];

  vecmat::vec3f Origin = vecmat::vec3f::vec(origin[0], origin[1], origin[2]);
  vecmat::vec3f Direction = vecmat::vec3f::vec(direction[0], direction[1], direction[2]);
  vecmat::vec3f V_0 = vecmat::vec3f::vec(v_0[0], v_0[1], v_0[2]);
  vecmat::vec3f V_1 = vecmat::vec3f::vec(v_1[0], v_1[1], v_1[2]);
  vecmat::vec3f V_2 = vecmat::vec3f::vec(v_2[0], v_2[1], v_2[2]);

  vecmat::mat33f A = vecmat::mat33f::mat(vecmat::vec3f::zero() - Direction, V_1 - V_0, V_2 - V_0);
  vecmat::vec3f b = vecmat::vec3f::vec(Origin - V_0);
  vecmat::mat33f A_0 = vecmat::mat33f::mat(Origin - V_0, V_1 - V_0, V_2 - V_0);
  vecmat::mat33f A_1 = vecmat::mat33f::mat(vecmat::vec3f::zero() - Direction, Origin - V_0, V_2 - V_0);
  vecmat::mat33f A_2 = vecmat::mat33f::mat(vecmat::vec3f::zero() - Direction, V_1 - V_0, Origin - V_0);

  float one_over_determinant = 1.0f / A.determinant();
  float t = A_0.determinant() * one_over_determinant;
  float uu = A_1.determinant() * one_over_determinant;
  float vv = A_2.determinant() * one_over_determinant;

  if (t >= ray.tNear && t <= ray.tFar && uu >= 0 && vv >= 0 && uu + vv <= 1) {
    ray.tFar = t;
    *primID = this->primID;
    *u = uu;
    *v = vv;
    return true;
  }

  return false;
}

void Triangle::fillIntersection(float distance, int primID, float u, float v,
  Intersection* intersection) const {
  // 该函数实际上不会被调用
  return;
}

//--- TriangleMesh ---
TriangleMesh::TriangleMesh(const Json& json) : Shape(json) {
  const auto& filepath = fetchRequired<std::string>(json, "file");
  meshData = MeshData::loadFromFile(filepath);
}

RTCGeometry TriangleMesh::getEmbreeGeometry(RTCDevice device) const {
  RTCGeometry geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

  float* vertexBuffer = (float*)rtcSetNewGeometryBuffer(
    geometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float),
    meshData->vertexCount);
  for (int i = 0; i < meshData->vertexCount; ++i) {
    Point3f vertex = transform.toWorld(meshData->vertexBuffer[i]);
    vertexBuffer[3 * i] = vertex[0];
    vertexBuffer[3 * i + 1] = vertex[1];
    vertexBuffer[3 * i + 2] = vertex[2];
  }

  unsigned* indexBuffer = (unsigned*)rtcSetNewGeometryBuffer(
    geometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
    3 * sizeof(unsigned), meshData->faceCount);
  for (int i = 0; i < meshData->faceCount; ++i) {
    indexBuffer[i * 3] = meshData->faceBuffer[i][0].vertexIndex;
    indexBuffer[i * 3 + 1] = meshData->faceBuffer[i][1].vertexIndex;
    indexBuffer[i * 3 + 2] = meshData->faceBuffer[i][2].vertexIndex;
  }
  rtcCommitGeometry(geometry);
  return geometry;
}

bool TriangleMesh::rayIntersectShape(Ray& ray, int* primID, float* u,
  float* v) const {
  //* 当使用embree加速时，该方法不会被调用
  int geomID = -1;
  return acceleration->rayIntersect(ray, &geomID, primID, u, v);
}

void TriangleMesh::fillIntersection(float distance, int primID, float u,
  float v, Intersection* intersection) const {
  //* todo 填充光线与三角网格求交得到的交点信息
  intersection->distance = distance;
  intersection->shape = this;
  //* 1. 在三角形内部用插值计算交点坐标
  //* 2. 在三角形内部用插值计算法线
  //* 3. 在三角形内部用插值计算纹理坐标
  //* 4. 在三角形内部用插值计算交点的切线和副切线

  Point3f v_0 = meshData->vertexBuffer[meshData->faceBuffer[primID][0].vertexIndex],
    v_1 = meshData->vertexBuffer[meshData->faceBuffer[primID][1].vertexIndex],
    v_2 = meshData->vertexBuffer[meshData->faceBuffer[primID][2].vertexIndex];
  intersection->position = Point3f((1 - u - v) * v_0[0] + u * v_1[0] + v * v_2[0],
    (1 - u - v) * v_0[1] + u * v_1[1] + v * v_2[1],
    (1 - u - v) * v_0[2] + u * v_1[2] + v * v_2[2]);

  Vector3f n_0 = meshData->normalBuffer[meshData->faceBuffer[primID][0].normalIndex],
    n_1 = meshData->normalBuffer[meshData->faceBuffer[primID][1].normalIndex],
    n_2 = meshData->normalBuffer[meshData->faceBuffer[primID][2].normalIndex];
  intersection->normal = normalize((1 - u - v) * n_0 + u * n_1 + v * n_2);

  Vector2f tex_0 = meshData->texcodBuffer[meshData->faceBuffer[primID][0].texcodIndex],
    tex_1 = meshData->texcodBuffer[meshData->faceBuffer[primID][1].texcodIndex],
    tex_2 = meshData->texcodBuffer[meshData->faceBuffer[primID][2].texcodIndex];
  intersection->texCoord = (1 - u - v) * tex_0 + u * tex_1 + v * tex_2;

  Vector3f tangent{ 1.f, 0.f, .0f };
  Vector3f bitangent;
  if (std::abs(dot(tangent, intersection->normal)) > .9f) {
    tangent = Vector3f(.0f, 1.f, .0f);
  }
  bitangent = normalize(cross(tangent, intersection->normal));
  tangent = normalize(cross(intersection->normal, bitangent));
  intersection->tangent = tangent;
  intersection->bitangent = bitangent;
}

void TriangleMesh::initInternalAcceleration() {
  acceleration = Acceleration::createAcceleration();
  int primCount = meshData->faceCount;
  for (int primID = 0; primID < primCount; ++primID) {
    int vtx0Idx = meshData->faceBuffer[primID][0].vertexIndex,
      vtx1Idx = meshData->faceBuffer[primID][1].vertexIndex,
      vtx2Idx = meshData->faceBuffer[primID][2].vertexIndex;
    std::shared_ptr<Triangle> triangle =
      std::make_shared<Triangle>(primID, vtx0Idx, vtx1Idx, vtx2Idx, this);
    acceleration->attachShape(triangle);
  }
  acceleration->build();
  // TriangleMesh的包围盒就是其内部加速结构的包围盒
  boundingBox = acceleration->boundingBox;
}
REGISTER_CLASS(TriangleMesh, "triangle")