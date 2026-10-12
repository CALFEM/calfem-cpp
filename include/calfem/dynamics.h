//
// CALFEM C++ library
// Copyright (c) 1999-2026 Division of Structural Mechanics, Lund University
//
// Written by Jonas Lindemann
//
// SPDX-License-Identifier: MIT
// See the LICENSE file in the project root for full license information.

// Time integration of first and second order systems.

#pragma once

#include <calfem/types.h>

#include <tuple>

namespace calfem {

// Solution history. Each column holds the values at one stored time.

struct TimeHistory {
    Matrix a;   // displacements (or the primary variable)
    Matrix da;  // first time derivative
    Matrix d2a; // second time derivative (step2 only)
};

// ---- step1: first order system C da + K a = f (e.g. transient heat flow) ----
//
//   modelhist, dofhist = step1(K, C, f, a0, bc, ip, times, dofs)
//   auto [modelhist, dofhist] = step1(K, C, f, a0, bc, ip, times, dofs);
//
// K, C   system matrices (ndof x ndof)
// f      load: ndof x 1 (constant) or ndof x (nstep+1) (one column per time step); empty or
//        all zero for no load
// a0     initial values (ndof)
// bc     prescribed values: [dof, value] rows (constant) or [dof, v0, v1, ..., vnstep]
//        (one column per time step); empty for none
// ip     [dt, T, alpha]: time step, total time and integration parameter (0 forward Euler,
//        0.5 trapezoidal rule (Crank-Nicolson), 1 backward Euler)
// times  times at which the whole solution is stored in modelhist; empty to store all steps
// dofs   dofs (1-based) whose history is stored in dofhist at every step; empty for none
//
// The number of steps is int(T/dt), or the number of time steps given in f or bc if smaller,
// as in calfem-python. modelhist stores a and da (one column per stored time), dofhist stores
// the history of the selected dofs (one row per dof, one column per step).
//
// Deviation from calfem-python: dofhist is recorded at every step also when times is empty
// (calfem-python then leaves all columns after the first at zero).
//
// Throws SolverError if C or the effective system matrix is singular.

std::tuple<TimeHistory, TimeHistory> step1(const Matrix &K, const Matrix &C, const Matrix &f, const ColVec &a0,
                                           const Matrix &bc, const RowVec &ip, const RowVec &times,
                                           const IntColVec &dofs);

// ---- step2: second order system M d2a + C da + K a = f (structural dynamics, Newmark) ----
//
//   auto [modelhist, dofhist] = step2(K, C, M, f, a0, da0, bc, ip, times, dofs);
//
// As step1, with M the mass matrix (C may be empty or zero for no damping), da0 the initial
// velocities and ip = [dt, T, alpha, delta] (Newmark parameters, e.g. alpha = 0.25,
// delta = 0.5 for the average acceleration method). modelhist and dofhist also store d2a.
// The accelerations of prescribed dofs are stored as zero, as in calfem-python.
//
// Throws SolverError if M or the effective system matrix is singular.

std::tuple<TimeHistory, TimeHistory> step2(const Matrix &K, const Matrix &C, const Matrix &M, const Matrix &f,
                                           const ColVec &a0, const ColVec &da0, const Matrix &bc, const RowVec &ip,
                                           const RowVec &times, const IntColVec &dofs);

// ---- gfunc: piecewise linear time function ----
//
//   ti, g = gfunc(G, dt)        auto [ti, g] = gfunc(G, dt);
//
// G = [[t0, g0], [t1, g1], ...] with increasing times. Returns the times ti = t0, t0 + dt, ...
// up to the last time in G and the linearly interpolated values g(ti).
//
// Deviation from calfem-python: the time points are t0 + i dt for all points not after the last
// time in G (np.arange may add a point beyond it because of rounding).

std::tuple<ColVec, ColVec> gfunc(const Matrix &G, double dt);

} // namespace calfem
