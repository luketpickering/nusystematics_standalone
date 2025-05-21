#include "nusystematics/systproviders/CVspread_tool.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "RwFramework/GSyst.h"

// GENIE-MC/Generator
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"

#include <cmath>

using namespace nusyst;
using namespace systtools;

CVspread::CVspread(fhicl::ParameterSet const &params)
    : IGENIESystProvider_tool(params) {}

std::vector<std::unique_ptr<TH3D>> geths(std::string const &loc) {
  TFile fin(loc.c_str(), "READ");
  std::vector<std::unique_ptr<TH3D>> hists;
  for (auto hn :
       {"CC0pi", "CC1pip", "CC1pim", "CC1pi0", "CC2cpi", "CCNpi", "CCOther",
        "NC0pi", "NC1pip", "NC1pim", "NC1pi0", "NC2cpi", "NCNpi", "NCOther"}) {
    hists.emplace_back(fin.Get<TH3D>(hn));
    hists.back()->SetDirectory(nullptr);
  }
  return hists;
}

void CVspread::LoadInputs() {
  auto inputs = tool_options.get<fhicl::ParameterSet>("inputs");

  ref_xs = geths(inputs.get<std::string>("reference"));

  for (auto &di : dial_infos) {
    di.alt_xs = geths(inputs.get<std::string>(di.prettyname));
  }
}

SystMetaData CVspread::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                         paramId_t firstId) {
  SystMetaData smd;
  SystParamHeader dial_variation_template;
  dial_variation_template.isSplineable = true;
  dial_variation_template.paramVariations = {0, 1};

  fhicl::ParameterSet inputs;
  inputs.put("ref", ps.get<std::string>("reference"));

  std::vector<std::string> alt_models;

  for (auto const &altm :
       ps.get<std::vector<fhicl::ParameterSet>>("alt_models")) {

    SystParamHeader phdr = dial_variation_template;
    phdr.prettyName = altm.get<std::string>("name");
    alt_models.push_back(phdr.prettyName);

    inputs.put(phdr.prettyName, altm.get<std::string>("input"));

    phdr.systParamId = firstId++;

    smd.push_back(phdr);
  }

  tool_options.put("inputs", inputs);
  tool_options.put("alt_models", alt_models);
  tool_options.put("verbosity_level", ps.get<int>("verbosity_level", 0));

  return smd;
}

bool CVspread::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {
  verbosity_level = tool_options.get<int>("verbosity_level", 0);

  std::vector<std::string> alt_models =
      tool_options.get<std::vector<std::string>>("alt_models");

  // grab the pre-parsed param headers object
  SystMetaData const &md = GetSystMetaData();

  for (auto altm : alt_models) {

    if (!HasParam(md, altm)) {
      if (verbosity_level > 1) {
        std::cout << "[INFO]: Don't have parameter " << altm
                  << " in SystMetaData. Skipping configuration." << std::endl;
      }
      continue;
    }

    auto pid = GetParamIndex(md, altm);

    if (verbosity_level > 1) {
      std::cout << "[INFO]: Have parameter " << altm
                << " in SystMetaData with ParamId: " << pid << ". Configuring."
                << std::endl;
    }
    dial_infos.push_back(DialInfo{altm, pid, {}});
  }
  LoadInputs();
  // returning cleanly
  return true;
}

namespace {
struct Topology {
  enum topo {
    kCC0pi = 0,
    kCC1pip,
    kCC1pim,
    kCC1pi0,
    kCC2cpi,
    kCCNpi,
    kCCOther,
    kNC0pi,
    kNC1pip,
    kNC1pim,
    kNC1pi0,
    kNC2cpi,
    kNCNpi,
    kNCOther,
    kInvalid
  };
};
} // namespace

Topology::topo get_reweight_topology(genie::EventRecord const &ev) {

  int nnu = 0;
  int nclep = 0;
  int nnuc = 0;
  int npip = 0;
  int npim = 0;
  int npi0 = 0;
  int ngamma = 0;
  int nother = 0;

  for (auto const &po : ev) {
    genie::GHepParticle const &p =
        dynamic_cast<genie::GHepParticle const &>(*po);

    if (p.Status() != genie::kIStStableFinalState) {
      continue;
    }
    switch (p.Pdg()) {
    case -12:
    case -14:
    case -16:
    case 12:
    case 14:
    case 16: {
      nnu++;
      break;
    }
    case -11:
    case -13:
    case -15:
    case 11:
    case 13:
    case 15: {
      nclep++;
      break;
    }
    case 2212:
    case 2112: {
      nnuc++;
      break;
    }
    case -211: {
      npim++;
      break;
    }
    case 211: {
      npip++;
      break;
    }
    case 111: {
      npi0++;
      break;
    }
    case 22: {
      ngamma++;
      break;
    }
    default: {
      if (p.Pdg() < 1000000000) {
        nother++;
      }
      break;
    }
    }
  }

  if (nclep > 1) {
    return Topology::kCCOther;
  }

  if (nnu > 1) {
    return Topology::kInvalid;
  }

  bool iscc = !((nnu == 1) && (nclep == 0));

  if (nother) {
    return iscc ? Topology::kCCOther : Topology::kNCOther;
  }

  if (ngamma) {
    return iscc ? Topology::kCCOther : Topology::kNCOther;
  }

  if (npi0) {
    return (npi0 > 1) ? (iscc ? Topology::kCCNpi : Topology::kNCNpi)
                      : (iscc ? Topology::kCC1pi0 : Topology::kNC1pi0);
  }

  if ((npip + npim)) {
    if ((npip + npim) > 2) {
      return iscc ? Topology::kCCNpi : Topology::kNCNpi;
    }
    if ((npip + npim) == 2) {
      return iscc ? Topology::kCC2cpi : Topology::kNC2cpi;
    }
    return (npip == 1) ? (iscc ? Topology::kCC1pip : Topology::kNC1pip)
                       : (iscc ? Topology::kCC1pim : Topology::kNC1pim);
  }

  return iscc ? Topology::kCC0pi : Topology::kNC0pi;
}

event_unit_response_t CVspread::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp;

  SystMetaData const &md = GetSystMetaData();

  genie::GHepParticle *ISLep = ev.Probe();

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();

  if (!FSLep || !ISLep) {
    return resp;
  }

  auto topo = get_reweight_topology(ev);

  if (topo == Topology::kInvalid) {
    return resp;
  }

  auto em_transfer_GeV = (*ISLep->P4() - *FSLep->P4());
  auto Q2_GeV = -em_transfer_GeV.Mag2();
  auto q_0_GeV = em_transfer_GeV.E();
  auto Enu_true = ISLep->P4()->E();

  constexpr static double m_p_GeV = 0.93827203;

  auto W_nuc_rest_GeV =
      std::sqrt((-Q2_GeV) + (2 * m_p_GeV * q_0_GeV) + (m_p_GeV * m_p_GeV));

  // get individual axis bins and check that the value is in range.
  int binX = ref_xs[topo]->GetXaxis()->FindBin(Enu_true);
  if ((binX == 0) || (binX == ref_xs[topo]->GetXaxis()->GetNbins() + 1)) {
    return resp;
  }
  int binY = ref_xs[topo]->GetYaxis()->FindBin(W_nuc_rest_GeV);
  if ((binY == 0) || (binY == ref_xs[topo]->GetYaxis()->GetNbins() + 1)) {
    return resp;
  }
  int binZ = ref_xs[topo]->GetZaxis()->FindBin(Q2_GeV);
  if ((binZ == 0) || (binZ == ref_xs[topo]->GetZaxis()->GetNbins() + 1)) {
    return resp;
  }

  int binGlobal = ref_xs[topo]->GetBin(binX, binY, binZ);

  // Get the bin content
  double ev_ref_xs = ref_xs[topo]->GetBinContent(binGlobal);

  // loop through and calculate weights
  for (auto const &di : dial_infos) {

    resp.push_back({di.pid, {}});

    double ev_xs_ratio = di.alt_xs[topo]->GetBinContent(binGlobal) / ev_ref_xs;

    resp.back().responses = {1.0, 1.0 + (ev_xs_ratio - 1.0)};
    if (verbosity_level > 3) {
      std::cout << "[DEBG]: For parameter " << di.prettyname << " at variation["
                << 1 << "] = " << 1
                << " calculated weight: " << resp.back().responses.back()
                << std::endl;
    }
  }

  return resp;
}
