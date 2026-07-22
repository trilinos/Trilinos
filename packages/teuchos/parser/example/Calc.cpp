// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <Teuchos_Language.hpp>
#include <Teuchos_Reader.hpp>
#include <Teuchos_MathExpr.hpp>

int main() {
  Teuchos::RCP<Teuchos::Reader> reader(Teuchos::MathExpr::new_calc_reader());
  for (std::string line; std::getline(std::cin, line);) {
    Teuchos::any result_any;
    try {
      reader->read_string(result_any, line, "input");
      double value = Teuchos::any_cast<double>(result_any);
      std::cout << value << '\n';
    } catch (const Teuchos::ParserFail& e) {
      std::cerr << e.what() << '\n';
    }
  }
}
