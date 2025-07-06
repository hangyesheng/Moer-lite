/**
 * @file  DisneyBSDF.h
 * @author Junping Yuan
 * @brief  Disney Principled BSDF.Highly inspired by https://cseweb.ucsd.edu/~tzli/cse272/homework1.pdf.
 * @version 0.1
 * @date 2022.10.13
 *
 * @copyright NJUMeta (c) 2022
 * www.njumeta.com
 *
 */
#include "FunctionLayer/Material/BxDF/BSDF.h"
#include "Material.h"
 // #include "BxDF/MicrofacetDistribution.h"
 // #include "FunctionLayer/Material/BxDF/DielectricBxDF.h"

struct DisneyDiffuse;
struct DisneyMetal;
struct DisneyClearCoat;
struct DisneyGlass;
struct DisneySheen;


class DisneyBSDF : public  BSDF {
public:
    DisneyBSDF(const Spectrum& baseColor, double specularTransmission, double metallic, double subsurface,
        double specular, double roughness, double specularTint, double anisotropic, double sheen,
        double sheenTint, double clearCoat, double clearCoatGloss, double eta);

protected:
    BSDFSampleResult sample(const Vector3f& out, const Vector2f& sample) const override;
    Spectrum f(const Vector3f& out, const Vector3f& in) const override;
    double  specularTransmission;
    double  metallic;
    double  specular;
    double  roughness;
    double  sheen;
    double  clearCoat;
    double  ior;

    std::unique_ptr<DisneyGlass> disneyGlass;
    std::unique_ptr<DisneyDiffuse> disneyDiffuse;
    std::unique_ptr<DisneyMetal> disneyMetal;
    std::unique_ptr<DisneyClearCoat> disneyClearCoat;
    std::unique_ptr<DisneySheen> disneySheen;

};

class DisneyMaterial : public Material {
public:
    DisneyMaterial(const Json& json);
    std::shared_ptr<BSDF> computeBSDF(const Intersection& intersection) const override;
protected:
    std::shared_ptr<Texture<Spectrum>> baseColor;
    std::shared_ptr<Texture<double>>  specularTransmission;
    std::shared_ptr<Texture<double>>  metallic;
    std::shared_ptr<Texture<double>>  subsurface;
    std::shared_ptr<Texture<double>>  specular;
    std::shared_ptr<Texture<double>>  roughness;
    std::shared_ptr<Texture<double>>  specularTint;
    std::shared_ptr<Texture<double>>  anisotropic;
    std::shared_ptr<Texture<double>>  sheen;
    std::shared_ptr<Texture<double>>  sheenTint;
    std::shared_ptr<Texture<double>>  clearCoat;
    std::shared_ptr<Texture<double>>  clearCoatGloss;

    double eta;
};