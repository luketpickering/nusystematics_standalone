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

std::map<int, std::set<int>> get_configurations(fhicl::ParameterSet const &ps) {
  std::map<int, std::set<int>> cfgs;
  for (auto const &cfg :
       ps.get<std::vector<fhicl::ParameterSet>>("configurations")) {
    cfgs[cfg.get<int>("tgt")].insert(cfg.get<int>("nu"));
  }
  return cfgs;
}

static std::vector<std::string> const Topologies = {
    "CC0pi", "CC1pip", "CC1pim", "CC1pi0", "CC2cpi", "CCGamma", "CCNpi", "CCOther",
    "NC0pi", "NC1pip", "NC1pim", "NC1pi0", "NC2cpi", "NCGamma", "NCNpi", "NCOther"};

std::map<size_t, std::map<size_t, std::vector<std::unique_ptr<TH3D>>>>
geths(std::string const &name, fhicl::ParameterSet const &ps, int verbosity_level) {

  auto input = ps.get<std::string>("input");
  auto cfgs = get_configurations(ps);

  if (verbosity_level > 1) {
    std::cout << "[INFO]: CVspread Reading histograms from file: " << input
              << std::endl;
  }

  TFile fin(input.c_str(), "READ");

  std::map<size_t, std::map<size_t, std::vector<std::unique_ptr<TH3D>>>> hists;

  for (auto const &[tgt, nupids] : cfgs) {
    std::map<size_t, std::vector<std::unique_ptr<TH3D>>> species;
    for (auto const &nupid : nupids) {
      std::vector<std::unique_ptr<TH3D>> topos;
      for (auto const &topo : Topologies) {
        auto hname = name + "_" + topo + "_" + std::to_string(nupid) + "_" +
                     std::to_string(tgt);

        if (verbosity_level > 1) {
          std::cout << "  - hist named: " << hname;
        }

        auto h = fin.Get<TH3D>(hname.c_str());
        if (h) {
          h->SetDirectory(nullptr);
          if (verbosity_level > 1) {
            std::cout << " exists!" << std::endl;
          }
        } else {
          if (verbosity_level > 1) {
            std::cout << " doesn't exist." << std::endl;
          }
        }
        topos.emplace_back(h);
      }
      species[nupid] = std::move(topos);
    }
    hists[tgt] = std::move(species);
  }
  return hists;
}

SystMetaData CVspread::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                         paramId_t firstId) {
  SystMetaData smd;
  SystParamHeader dial_variation_template;
  dial_variation_template.isSplineable = true;
  dial_variation_template.paramVariations = {0, 1};

  fhicl::ParameterSet inputs;

  auto ref_obj = ps.get<fhicl::ParameterSet>("reference");
  inputs.put("reference", ref_obj);

  auto ref_cfgs = get_configurations(ref_obj);

  std::vector<fhicl::ParameterSet> alternate_models;

  for (auto altm :
       ps.get<std::vector<fhicl::ParameterSet>>("alternate_models")) {

    auto name = altm.get<std::string>("name");

    std::map<int, std::set<int>> alt_cfgs;

    if (altm.has_key("configurations")) {
      alt_cfgs = get_configurations(altm);
      for (auto const &[tgt, nupids] : alt_cfgs) {
        if (!ref_cfgs.count(tgt)) {
          std::stringstream ss;
          ss << "Alternate model: " << name << " provides target: " << tgt
             << " for which we have no reference xsec.";
          throw std::runtime_error(ss.str());
        }
        for (auto const &nupid : nupids) {
          if (!ref_cfgs[tgt].count(nupid)) {
            std::stringstream ss;
            ss << "Alternate model: " << name
               << " provides neutrino species: " << nupid
               << " on target: " << tgt
               << " for which we have no reference xsec.";
            throw std::runtime_error(ss.str());
          }
        }
      }
    } else {
      alt_cfgs = ref_cfgs;
      altm.put("configurations",
               ref_obj.get<std::vector<fhicl::ParameterSet>>("configurations"));
    }

    for (auto const &[tgt, nupids] : alt_cfgs) {
      for (auto const &nupid : nupids) {
        for (auto const &topo : Topologies) {
          SystParamHeader phdr = dial_variation_template;
          phdr.prettyName = name + "_" + topo + "_" + std::to_string(nupid) +
                            "_" + std::to_string(tgt);
          phdr.systParamId = firstId++;
          smd.push_back(phdr);
        }
      }
    }

    alternate_models.push_back(altm);
  }

  inputs.put("alternate_models", alternate_models);

  tool_options.put("inputs", inputs);
  tool_options.put("verbosity_level", ps.get<int>("verbosity_level", 0));

  return smd;
}

bool CVspread::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {
  verbosity_level = tool_options.get<int>("verbosity_level", 0);

  auto inputs = tool_options.get<fhicl::ParameterSet>("inputs");

  auto ref_ps = inputs.get<fhicl::ParameterSet>("reference");

  auto name = ref_ps.get<std::string>("name");
  auto ref_xs = geths(name, ref_ps, verbosity_level);

  // grab the pre-parsed param headers object
  SystMetaData const &md = GetSystMetaData();

  for (auto altm :
       inputs.get<std::vector<fhicl::ParameterSet>>("alternate_models")) {

    auto name = altm.get<std::string>("name");
    auto altm_xs = geths(name, altm, verbosity_level);

    for (auto &[tgt, nupids] : altm_xs) {
      for (auto &[nupid, hists] : nupids) {
        for (size_t h_it = 0; h_it < hists.size(); ++h_it) {

          if (!hists[h_it]) {
            continue;
          }

          auto topo = Topologies[h_it];

          auto dial_prettyname = name + "_" + topo + "_" +
                                 std::to_string(nupid) + "_" +
                                 std::to_string(tgt);

          if (!HasParam(md, dial_prettyname)) {
            if (verbosity_level > 1) {
              std::cout << "[INFO]: Don't have parameter " << dial_prettyname
                        << " in SystMetaData. Skipping configuration."
                        << std::endl;
            }
            continue;
          }

          auto pid = GetParamIndex(md, dial_prettyname);

          if (verbosity_level > 1) {
            std::cout << "[INFO]: Have parameter " << dial_prettyname
                      << " in SystMetaData with ParamId: " << pid
                      << ". Configuring." << std::endl;
          }

          dial_infos.emplace_back(dial_prettyname, pid, tgt, nupid, h_it,
                                  std::move(hists[h_it]));
          auto ref = ref_xs[tgt][nupid][h_it].get();
          dial_infos.back().weights->Divide(ref);
        }
      }
    }
  }

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
    kCCGamma,
    kCCOther,
    kNC0pi,
    kNC1pip,
    kNC1pim,
    kNC1pi0,
    kNC2cpi,
    kNCNpi,
    kNCGamma,
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
      if(p.P4()->E() > 0.01){
        ngamma++;
      }
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
    return iscc ? Topology::kCCGamma : Topology::kNCGamma;
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

  if (!FSLep || !ISLep || !ev.TargetNucleus()) {
    return resp;
  }

  auto TargetPDG = ev.TargetNucleus()->Pdg();

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

  // loop through and calculate weights
  for (auto const &di : dial_infos) {

    if ((di.target_pid != TargetPDG) || (di.nu_pid != ISLep->Pdg()) ||
        (di.topology != topo) || (!di.weights)) {
      continue;
    }

    resp.push_back({di.param_id, {1.0, 1.0}});

    // get individual axis bins and check that the value is in range.
    int binX = di.weights->GetXaxis()->FindFixBin(Enu_true);
    if ((binX == 0) || (binX == di.weights->GetXaxis()->GetNbins() + 1)) {
      continue;
    }
    int binY = di.weights->GetYaxis()->FindFixBin(W_nuc_rest_GeV);
    if ((binY == 0) || (binY == di.weights->GetYaxis()->GetNbins() + 1)) {
      continue;
    }
    int binZ = di.weights->GetZaxis()->FindFixBin(Q2_GeV);
    if ((binZ == 0) || (binZ == di.weights->GetZaxis()->GetNbins() + 1)) {
      continue;
    }

    auto w = di.weights->GetBinContent(binX, binY, binZ);

    resp.back().responses[1] = std::min(std::max(w, 0.0), 10.0);

    if (verbosity_level > 3) {
      std::cout << "[DEBG]: For parameter " << di.prettyname << " at variation["
                << 1 << "] = " << 1
                << " calculated weight: " << resp.back().responses.back()
                << std::endl;
    }
  }

  return resp;
}
