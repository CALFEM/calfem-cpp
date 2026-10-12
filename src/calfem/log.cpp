//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

#include <calfem/log.h>

#include "detail/log.h"

namespace calfem {

namespace {

LogFunction &logFunction()
{
    static LogFunction func;
    return func;
}

} // namespace

void setLogFunction(LogFunction func)
{
    logFunction() = std::move(func);
}

namespace detail {

void logMessage(const std::string &message)
{
    auto &func = logFunction();

    if (func)
        func(message);
}

std::string &lastError()
{
    thread_local std::string error;
    return error;
}

void logError(const std::string &message)
{
    lastError() = message;
    logMessage(message);
}

} // namespace detail

} // namespace calfem
