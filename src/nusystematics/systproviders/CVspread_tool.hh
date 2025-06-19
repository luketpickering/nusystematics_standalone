#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "TFile.h"
#include "TH3D.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include <array>
#include <memory>
#include <string>

class CVspread : public nusyst::IGENIESystProvider_tool {

public:
  explicit CVspread(fhicl::ParameterSet const &);

  //'First' configuration step: tool configuration
  // - takes 'arbitrary' FHiCL configuration and
  //   produces SystMetaData object which can later be used to configure a
  //   specific set of parameter values to be calculated
  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  //'Second' configuration step: parameter headers
  // - Reads the preconstructed SystMetaData produced by BuildSystMetaData
  //   to configure an instance of this class to calculate weights
  // - Recieves a copy of the tool_options instance constructed by
  //   BuildSystMetaData as an argument
  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  // Used to pass arbitrary FHiCL options from the tool configuration to the
  //   parameter headers.
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  // Parameter-specific implementation goes in here
  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  // Can add as much or as little stateful information here for use when
  // representing this instance as a string.
  std::string AsString() { return "CVspread"; }

  ~CVspread() {}

private:
  // arbitrary additional configuration from the tool configuration/parameter
  // headers can be storeds here
  fhicl::ParameterSet tool_options;

  struct DialInfo {

    DialInfo(std::string prettyname_, systtools::paramId_t param_id_,
             size_t target_pid_, size_t nu_pid_, size_t topology_,
             std::unique_ptr<TH3D> &&weights_)
        : prettyname(prettyname_), param_id(param_id_), target_pid(target_pid_),
          nu_pid(nu_pid_), topology(topology_), weights(std::move(weights_)) {}
    std::string prettyname;
    systtools::paramId_t param_id;
    size_t target_pid;
    size_t nu_pid;
    size_t topology;
    std::unique_ptr<TH3D> weights;
  };
  std::vector<DialInfo> dial_infos;

  // configurable verbosity as an example of some arbitrary systprovider
  // configuration
  int verbosity_level;
};