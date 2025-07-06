#include "MicrofacetDistribution.h"
// #include "CoreLayer/Geometry/Frame.h"

#include "BSDF.h"
// #include "FastMath.h"
#include <ResourceLayer/JsonUtil.h>

std::shared_ptr<MicrofacetDistribution> LoadDistributionFromJson(const Json& json) {
    std::string distribStr = json["distribution"];
    if (distribStr == "ggx")
        return std::make_shared<GGXDistribution>(true);
    {}
}

MicrofacetDistribution::~MicrofacetDistribution() noexcept {}

double MicrofacetDistribution::Pdf(const Vector3f& wo, const Vector3f& wh, const Vector2f& alphaXY) const {
    if (sampleVisibleArea)
        return D(wh, alphaXY) * G1(wo, alphaXY) * absDot(wo, wh) / AbsCosTheta(wo);
    else
        return D(wh, alphaXY) * AbsCosTheta(wh);
}

double GGXDistribution::G(const Vector3f& wo, const Vector3f& wi, const Vector2f& alphaXY) const {
    return G1(wo, alphaXY) * G1(wi, alphaXY);
}

double GGXDistribution::roughnessToAlpha(double roughness) const {
    return roughness;
}

double GGXDistribution::D(const Vector3f& wh, const Vector2f& alphaXY) const {
    double alphaX = alphaXY.x();
    double alphaY = alphaXY.y();
    double ax2 = alphaX * alphaX;
    double ay2 = alphaY * alphaY;
    Vector3f wh2 = wh * wh;
    double D = M_PI * alphaX * alphaY * pow(wh2.x() / ax2 + wh2.y() / ay2 + wh2.z(), 2);
    return 1 / D;
}

double GGXDistribution::Lambda(const Vector3f& w, const Vector2f& alphaXY) const {
    double ax2 = alphaXY.x() * alphaXY.x();
    double ay2 = alphaXY.y() * alphaXY.y();
    Vector3f v2 = w * w;
    double Lambda = (-1 + sqrt(1 + (v2.x() * ax2 + v2.y() * ay2) / v2.z())) / 2;
    return Lambda;
}

Vector3f GGXDistribution::Sample_wh(const Vector3f& wo, Vector2f u, const Vector2f& alphaXY) const {
    u = Vector2f(
        std::max(0.01f, std::min(0.99f, u.x())),
        std::max(0.01f, std::min(0.99f, u.y()))
    );
    double alphaX = alphaXY.x(), alphaY = alphaXY.y();
    if (CosTheta(wo) < 0) {
        return Sample_wh(-wo, u, alphaXY);
    }
    if (sampleVisibleArea) {
        // see https://jcgt.org/published/0007/04/01/slides.pdf
        //  Transform the incoming direction to the "hemisphere configuration".
        Vector3f hemisphereDirOut = normalize(Vector3f(alphaX * wo.x(), alphaY * wo.y(), wo.z()));
        // Parameterization of the projected area of a hemisphere.
        double r = sqrt(u.x());
        double phi = 2 * M_PI * u.y();
        float t1 = r * cos(phi);
        float t2 = r * sin(phi);
        // Vertically scale the position of a sample to account for the projection.
        double s = (1 + hemisphereDirOut.z()) / 2;
        t2 = (1 - s) * sqrt(1 - t1 * t1) + s * t2;
        // Point in the disk space
        Vector3f diskN{ t1, t2, sqrt(std::max(0.0f, 1 - t1 * t1 - t2 * t2)) };
        // Reprojection onto hemisphere -- we get our sampled normal in hemisphere space.
        Vector3f T1 = normalize(Vector3f(-hemisphereDirOut.y(), hemisphereDirOut.x(), 0));
        Vector3f T2 = cross(hemisphereDirOut, T1);
        Vector3f hemisphereN = t1 * T1 + t2 * T2 + diskN.z() * hemisphereDirOut;

        // Transforming the normal back to the ellipsoid configuration
        return normalize(Vector3f(alphaX * hemisphereN.x(), alphaY * hemisphereN.y(), std::max(0.0f, hemisphereN.z())));
    }
    else {
        double cosTheta, phi = (2 * M_PI) * u[1];
        if (alphaX == alphaY) {
            double tanTheta2 = alphaX * alphaY * u[0] / (1.0f - u[0]);
            cosTheta = 1 / fm::sqrt(1 + tanTheta2);
        }
        else {
            phi =
                fm::atan(alphaY / alphaX * fm::tan(2 * M_PI * u[1] + .5f * M_PI));
            if (u[1] > .5f) phi += M_PI;
            double sinPhi = fm::sin(phi), cosPhi = fm::cos(phi);
            const double alphaX2 = alphaX * alphaX, alphaY2 = alphaY * alphaY;
            const double alpha2 =
                1 / (cosPhi * cosPhi / alphaX2 + sinPhi * sinPhi / alphaY2);
            double tanTheta2 = alpha2 * u[0] / (1 - u[0]);
            cosTheta = 1 / fm::sqrt(1 + tanTheta2);
        }
        double sinTheta =
            fm::sqrt(std::max((double)0., (double)1. - cosTheta * cosTheta));
        Vector3f wh = Vector3f(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
        return wh;
    }
}

std::string GGXDistribution::ToString() const {
    return std::string();
}
