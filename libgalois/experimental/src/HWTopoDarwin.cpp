/*
 * This file belongs to the Galois project, a C++ library for exploiting parallelism.
 * The code is being released under the terms of the 3-Clause BSD License (a
 * copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

#include "galois/substrate/HWTopo.h"
#include "galois/gIO.h"
#include <sys/types.h>
#include <sys/sysctl.h>

using namespace galois::substrate;

namespace {

struct Policy {
  // number of "real" processors
  uint32_t numCpus;

  Policy() {
    size_t bufSize = sizeof(numCpus);
    if (sysctlbyname("hw.activecpu", &numCpus, &bufSize, NULL, 0) == -1) {
      GALOIS_SYS_DIE("Error querying number of cpus");
    } else if (bufSize != sizeof(numCpus)) {
      GALOIS_SYS_DIE("Error querying number of cpus");
    }
  }
};

static Policy& getPolicy() {
  static Policy A;
  return A;
}

} // namespace


std::pair<machineTopoInfo, std::vector<threadTopoInfo>>
galois::substrate::getHWTopo() {
  machineTopoInfo retMTI;

  retMTI.maxSockets   = getPolicy().numCpus * 2;
  retMTI.maxThreads   = getPolicy().numCpus * 2;
  retMTI.maxCores     = getPolicy().numCpus;
  retMTI.maxNumaNodes = 1;

  std::vector<threadTopoInfo> retTTI;
  retTTI.reserve(retMTI.maxThreads);

   for (unsigned i = 0; i < retMTI.maxThreads; ++i) {
     unsigned pid = i;
     retTTI.push_back(threadTopoInfo{i, pid, pid, pid, pid});
   }
  return std::make_pair(retMTI, retTTI);
}


//! binds current thread to OS HW context "proc"
bool galois::substrate::bindThreadSelf(unsigned osContext) {
  galois::gWarn("No cpu affinity on Cygwin.  Performance will be bad.");
  return false;
}
