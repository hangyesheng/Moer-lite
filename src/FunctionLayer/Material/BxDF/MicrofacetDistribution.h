/**
 * @file MicrofacetDistribution.h
 * @author JunPing.Yuan
 * @brief some MicrofactDistirbution models
 * @version 0.1
 * @date 2022-10-5
 *
 * @copyright NJUMeta (c) 2022
 * www.njumeta.com
 *
 */
#pragma once

#include "CoreLayer/Math/Geometry.h"
 // #include "CoreLayer/Adapter/JsonUtil.h"

class MicrofacetDistribution {
public:
    // MicrofacetDistribution Public Methods
    virtual ~MicrofacetDistribution();
    // https://www.pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models
    virtual double roughnessToAlpha(double roughness) const = 0;
    virtual double D(const Vector3f& wh, const Vector2f& alphaXY) const = 0;
    virtual double Lambda(const Vector3f& w, const Vector2f& alphaXY) const = 0;
    double G1(const Vector3f& w, const Vector2f& alphaXY) const {
        //    if (Dot(w, wh) * CosTheta(w) < 0.) return 0.;
        return 1 / (1 + Lambda(w, alphaXY));
    }
    virtual double G(const Vector3f& wo, const Vector3f& wi, const Vector2f& alphaXY) const {
        return G1(wo, alphaXY) * G1(wi, alphaXY);
        return 1 / (1 + Lambda(wo, alphaXY) + Lambda(wi, alphaXY));
    }
    virtual Vector3f Sample_wh(const Vector3f& wo, Vector2f u, const Vector2f& alphaXY) const = 0;
    double Pdf(const Vector3f& wo, const Vector3f& wh, const Vector2f& alphaXY) const;
    virtual std::string ToString() const = 0;
    inline bool sampleVisible() {
        return sampleVisibleArea;
    }

protected:
    // MicrofacetDistribution Protected Methods
    MicrofacetDistribution(bool sampleVisibleArea)
        : sampleVisibleArea(sampleVisibleArea) {
    }

    // MicrofacetDistribution Protected Data
    const bool sampleVisibleArea;
};


class GGXDistribution : public MicrofacetDistribution {
public:
    GGXDistribution(bool sampleVis = false) : MicrofacetDistribution(sampleVis) {}
    double roughnessToAlpha(double roughness) const override;

    double D(const Vector3f& wh, const Vector2f& alphaXY) const override;

    double G(const Vector3f& wo, const Vector3f& wi, const Vector2f& alphaXY) const override;

    double Lambda(const Vector3f& w, const Vector2f& alphaXY) const override;

    Vector3f Sample_wh(const Vector3f& wo, Vector2f u, const Vector2f& alphaXY) const override;

    std::string ToString() const override;
};

