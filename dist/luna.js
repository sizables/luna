(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
module.exports = require('./lib/lune.js');

},{"./lib/lune.js":2}],2:[function(require,module,exports){
/**
 * This library calculates the current phase of the moon
 * as well as finds the dates of the recent moon phases.
 *
 * Ported from python version found here:
 * https://bazaar.launchpad.net/~keturn/py-moon-phase/trunk/annotate/head:/moon.py
 *
 * Author: Ryan Seys (https://github.com/ryanseys)
 */

(function() {
  // Phases of the moon & precision

  var PRECISION = 0.05;
  var NEW = 0 / 4.0;
  var FIRST = 1 / 4.0;
  var FULL = 2 / 4.0;
  var LAST = 3 / 4.0;
  var NEXTNEW = 4 / 4.0;

  /**
   * Gets the Julian value from a date object.
   * Source: http://javascript.about.com/library/bljulday.htm
   * @return {Number} Julian number representation of the date.
   */
  Date.prototype.getJulian = function() {
    return (this.valueOf() / 86400000) - (this.getTimezoneOffset() / 1440) + 2440587.5;
  };

  /**
   * Converts a Number in Julian date form to a Date object.
   * Source: http://blog.bahrenburgs.com/2011/01/javascript-julian-day-conversions.html
   */
  Number.prototype.Julian2Date = function(inUTC) {
    var X = parseFloat(this)+0.5;
    var Z = Math.floor(X); //Get day without time
    var F = X - Z; //Get time
    var Y = Math.floor((Z-1867216.25)/36524.25);
    var A = Z+1+Y-Math.floor(Y/4);
    var B = A+1524;
    var C = Math.floor((B-122.1)/365.25);
    var D = Math.floor(365.25*C);
    var G = Math.floor((B-D)/30.6001);
    //must get number less than or equal to 12)
    var month = (G<13.5) ? (G-1) : (G-13);
    //if Month is January or February, or the rest of year
    var year = (month<2.5) ? (C-4715) : (C-4716);
    month -= 1; //Handle JavaScript month format
    var UT = B-D-Math.floor(30.6001*G)+F;
    var day = Math.floor(UT);
    //Determine time
    UT -= Math.floor(UT);
    UT *= 24;
    var hour = Math.floor(UT);
    UT -= Math.floor(UT);
    UT *= 60;
    var minute = Math.floor(UT);
    UT -= Math.floor(UT);
    UT *= 60;
    var second = Math.round(UT);

    if (inUTC) {
      return new Date(Date.UTC(year, month, day, hour, minute, second));
    } else {
      return new Date(year, month, day, hour, minute, second);
    }
  };

  /**
   * Astronomical Constants
   * @type {Object}
   */
  const c = {
    // JDN stands for Julian Day Number
    // Angles here are in degrees

    // 1980 January 0.0 in JDN
    // XXX: DateTime(1980).jdn yields 2444239.5 -- which one is right?
    epoch: 2444238.5,

    // Ecliptic longitude of the Sun at epoch 1980.0
    ecliptic_longitude_epoch: 278.833540,

    // Ecliptic longitude of the Sun at perigee
    ecliptic_longitude_perigee: 282.596403,

    // Eccentricity of Earth's orbit
    eccentricity: 0.016718,

    // Semi-major axis of Earth's orbit, in kilometers
    sun_smaxis: 1.49585e8,

    // Sun's angular size, in degrees, at semi-major axis distance
    sun_angular_size_smaxis: 0.533128,

    // Elements of the Moon's orbit, epoch 1980.0

    // Moon's mean longitude at the epoch
    moon_mean_longitude_epoch: 64.975464,
    // Mean longitude of the perigee at the epoch
    moon_mean_perigee_epoch: 349.383063,

    // Mean longitude of the node at the epoch
    node_mean_longitude_epoch: 151.950429,

    // Inclination of the Moon's orbit
    moon_inclination: 5.145396,

    // Eccentricity of the Moon's orbit
    moon_eccentricity: 0.054900,

    // Moon's angular size at distance a from Earth
    moon_angular_size: 0.5181,

    // Semi-mojor axis of the Moon's orbit, in kilometers
    moon_smaxis: 384401.0,
    // Parallax at a distance a from Earth
    moon_parallax: 0.9507,

    // Synodic month (new Moon to new Moon), in days
    synodic_month: 29.53058868,

    // Base date for E. W. Brown's numbered series of lunations (1923 January 16)
    lunations_base: 2423436.0,

    // #Properties of the Earth
    earth_radius: 6378.16
  };

  function fixangle(a) {
    return a - 360.0 * Math.floor(a/360.0);
  }

  /**
   * Convert degrees to radians
   * @param  {Number} d Angle in degrees
   * @return {Number}   Angle in radians
   */
  function torad(d) {
    return d * Math.PI / 180.0;
  }

  /**
   * Convert radians to degrees
   * @param  {Number} r Angle in radians
   * @return {Number}   Angle in degrees
   */
  function todeg(r) {
    return r * 180.0 / Math.PI;
  }

  function dsin(d) {
    return Math.sin(torad(d));
  }

  function dcos(d) {
    return Math.cos(torad(d));
  }

  /**
   * Solve the equation of Kepler.
   */
  function kepler(m, ecc) {
    var epsilon = 1e-6;

    m = torad(m);
    var e = m;
    while(1) {
      var delta = e - ecc * Math.sin(e) - m;
      e = e - delta / (1.0 - ecc * Math.cos(e));

      if (Math.abs(delta) <= epsilon) {
        break;
      }
    }

    return e;
  }

  /**
   * Finds the phase information for specific date.
   * @param  {Date} phase_date Date to get phase information of.
   * @return {Object}          Phase data
   */
  function phase(phase_date) {
    if(!phase_date) {
      phase_date = (new Date()).getJulian();
    }
    else {
      phase_date = phase_date.getJulian();
    }

    var day = phase_date - c.epoch;

    // Mean anomaly of the Sun
    var N = fixangle((360/365.2422) * day);
    //Convert from perigee coordinates to epoch 1980
    var M = fixangle(N + c.ecliptic_longitude_epoch - c.ecliptic_longitude_perigee);

    // Solve Kepler's equation
    var Ec = kepler(M, c.eccentricity);
    Ec = Math.sqrt((1 + c.eccentricity) / (1 - c.eccentricity)) * Math.tan(Ec/2.0);
    // True anomaly
    Ec = 2 * todeg(Math.atan(Ec));
    // Suns's geometric ecliptic longuitude
    var lambda_sun = fixangle(Ec + c.ecliptic_longitude_perigee);

    // Orbital distance factor
    var F = ((1 + c.eccentricity * Math.cos(torad(Ec))) / (1 - Math.pow(c.eccentricity, 2)));

    // Distance to Sun in km
    var sun_dist = c.sun_smaxis / F;
    var sun_angular_diameter = F * c.sun_angular_size_smaxis;

    // Calculation of the Moon's position

    // Moon's mean longitude
    var moon_longitude = fixangle(13.1763966 * day + c.moon_mean_longitude_epoch);

    // Moon's mean anomaly
    var MM = fixangle(moon_longitude - 0.1114041 * day - c.moon_mean_perigee_epoch);

    // Moon's ascending node mean longitude
    // MN = fixangle(c.node_mean_longitude_epoch - 0.0529539 * day)

    var evection = 1.2739 * Math.sin(torad(2*(moon_longitude - lambda_sun) - MM));

    // Annual equation
    var annual_eq = 0.1858 * Math.sin(torad(M));

    // Correction term
    var A3 = 0.37 * Math.sin(torad(M));

    var MmP = MM + evection - annual_eq - A3;

    // Correction for the equation of the centre
    var mEc = 6.2886 * Math.sin(torad(MmP));

    // Another correction term
    var A4 = 0.214 * Math.sin(torad(2 * MmP));

    // Corrected longitude
    var lP = moon_longitude + evection + mEc - annual_eq + A4;

    // Variation
    var variation = 0.6583 * Math.sin(torad(2*(lP - lambda_sun)));

    // True longitude
    var lPP = lP + variation;

    // Calculation of the phase of the Moon

    // Age of the Moon, in degrees
    var moon_age = lPP - lambda_sun;

    // Phase of the Moon
    var moon_phase = (1 - Math.cos(torad(moon_age))) / 2.0;

    // Calculate distance of Moon from the centre of the Earth
    var moon_dist = (c.moon_smaxis * (1 - Math.pow(c.moon_eccentricity,2))) / (1 + c.moon_eccentricity * Math.cos(torad(MmP + mEc)));

    // Calculate Moon's angular diameter
    var moon_diam_frac = moon_dist / c.moon_smaxis;
    var moon_angular_diameter = c.moon_angular_size / moon_diam_frac;

    // Calculate Moon's parallax (unused?)
    // moon_parallax = c.moon_parallax / moon_diam_frac

    var res = {
      'phase': fixangle(moon_age) / 360.0,
      'illuminated': moon_phase,
      'age': c.synodic_month * fixangle(moon_age) / 360.0,
      'distance': moon_dist,
      'angular_diameter': moon_angular_diameter,
      'sun_distance': sun_dist,
      'sun_angular_diameter': sun_angular_diameter
    };

    return res;
  }

  /**
   * Find time of phases of the moon which surround the current date.
   * Five phases are found, starting and ending with the new moons
   * which bound the current lunation.
   * @param  {Date} sdate Date to start hunting from (defaults to current date)
   * @return {Object}     Object containing recent past and future phases
   */
  function phase_hunt(sdate) {
    if(!sdate) {
      sdate = new Date();
    }

    var adate = new Date(sdate.valueOf()); // today!
    var x = 45; // go back 45 days!
    adate.setDate(adate.getDate() - x);

    var k1 = Math.floor((adate.getFullYear() + ((adate.getMonth()) * (1.0/12.0)) - 1900) * 12.3685);
    var nt1 = meanphase(adate, k1);
    adate = nt1;

    sdate = sdate.getJulian();
    var k2;
    while(1) {
      adate = adate + c.synodic_month;
      k2 = k1 + 1;
      var nt2 = meanphase(adate, k2);
      if(nt1 <= sdate && sdate < nt2) {
        break;
      }
      nt1 = nt2;
      k1 = k2;
    }
    var ks = [k1, k1, k1, k1, k2];
    var tphases = [NEW, FIRST, FULL, LAST, NEW];
    var phase_names = ['new_date', 'q1_date', 'full_date', 'q3_date', 'nextnew_date'];
    var phases = {};

    for (var i = 0; i < ks.length; i++) {
      phases[phase_names[i]] = truephase(ks[i], tphases[i]);
    }

    return phases;
  }

  /**
   * Given a K value used to determine the mean phase of the new
   * moon, and a phase selector (0.0, 0.25, 0.5, 0.75), obtain the
   * true, corrected phase time.
   * @param  {[type]} k      [description]
   * @param  {[type]} tphase [description]
   * @return {[type]}        [description]
   */
  function truephase(k, tphase) {

    var apcor = false;

    // add phase to new moon time
    k = k + tphase;
    // Time in Julian centuries from 1900 January 0.5
    var t = k / 1236.85;

    var t2 = t * t;
    var t3 = t2 * t;

    // Mean time of phase
    var pt = (
      2415020.75933 + c.synodic_month * k + 0.0001178 * t2 -
      0.000000155 * t3 + 0.00033 * dsin(166.56 + 132.87 * t -
      0.009173 * t2)
    );

    // Sun's mean anomaly
    var m = 359.2242 + 29.10535608 * k - 0.0000333 * t2 - 0.00000347 * t3;

    // Moon's mean anomaly
    var mprime = 306.0253 + 385.81691806 * k + 0.0107306 * t2 + 0.00001236 * t3;

    // Moon's argument of latitude
    var f = 21.2964 + 390.67050646 * k - 0.0016528 * t2 - 0.00000239 * t3;

    if ((tphase < 0.01) || (Math.abs(tphase - 0.5) < 0.01)) {

      // Corrections for New and Full Moon
      pt = pt + (
        (0.1734 - 0.000393 * t) * dsin(m) +
        0.0021 * dsin(2 * m) -
        0.4068 * dsin(mprime) +
        0.0161 * dsin(2 * mprime) -
        0.0004 * dsin(3 * mprime) +
        0.0104 * dsin(2 * f) -
        0.0051 * dsin(m + mprime) -
        0.0074 * dsin(m - mprime) +
        0.0004 * dsin(2 * f + m) -
        0.0004 * dsin(2 * f - m) -
        0.0006 * dsin(2 * f + mprime) +
        0.0010 * dsin(2 * f - mprime) +
        0.0005 * dsin(m + 2 * mprime)
      );

      apcor = true;
    }
    else if ((Math.abs(tphase - 0.25) < 0.01) || (Math.abs(tphase - 0.75) < 0.01)) {
        pt = pt + (
          (0.1721 - 0.0004 * t) * dsin(m) +
          0.0021 * dsin(2 * m) -
          0.6280 * dsin(mprime) +
          0.0089 * dsin(2 * mprime) -
          0.0004 * dsin(3 * mprime) +
          0.0079 * dsin(2 * f) -
          0.0119 * dsin(m + mprime) -
          0.0047 * dsin(m - mprime) +
          0.0003 * dsin(2 * f + m) -
          0.0004 * dsin(2 * f - m) -
          0.0006 * dsin(2 * f + mprime) +
          0.0021 * dsin(2 * f - mprime) +
          0.0003 * dsin(m + 2 * mprime) +
          0.0004 * dsin(m - 2 * mprime) -
          0.0003 * dsin(2 * m + mprime)
        );
      if (tphase < 0.5) {
          //  First quarter correction
          pt = pt + 0.0028 - 0.0004 * dcos(m) + 0.0003 * dcos(mprime);
      }
      else {
          //  Last quarter correction
          pt = pt + -0.0028 + 0.0004 * dcos(m) - 0.0003 * dcos(mprime);
      }
      apcor = true;
    }

    if (!apcor) {
      console.log("TRUEPHASE called with invalid phase selector ", tphase);
    }

    return pt.Julian2Date(true);
  }

  /**
   * Calculates time of the mean new Moon for a given base date.
   * This argument K to this function is the precomputed synodic month
   * index, given by:
   *   K = (year - 1900) * 12.3685
   * where year is expressed as a year and fractional year.
   * @param  {Date} sdate   Start date
   * @param  {[type]} k     [description]
   * @return {[type]}       [description]
   */
  function meanphase(sdate, k) {

    // Time in Julian centuries from 1900 January 12 noon
    var delta_t = (sdate - (new Date(1900,0,1,12))) / (1000*60*60*24);
    var t = delta_t / 36525;

    // square for frequent use
    var t2 = t * t;
    // and cube
    var t3 = t2 * t;

    nt1 = (
      2415020.75933 + c.synodic_month * k + 0.0001178 * t2 -
      0.000000155 * t3 + 0.00033 * dsin(166.56 + 132.87 * t -
      0.009173 * t2)
    );

    return nt1;
  }

  module.exports = {
   'phase_hunt': phase_hunt,
   'phase': phase
  };
})();

},{}],3:[function(require,module,exports){
// (c) Dean McNamee <dean@gmail.com>, 2013.
//
// https://github.com/deanm/omggif
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// omggif is a JavaScript implementation of a GIF 89a encoder and decoder,
// including animation and compression.  It does not rely on any specific
// underlying system, so should run in the browser, Node, or Plask.

function GifWriter(buf, width, height, gopts) {
  var p = 0;

  var gopts = gopts === undefined ? { } : gopts;
  var loop_count = gopts.loop === undefined ? null : gopts.loop;
  var global_palette = gopts.palette === undefined ? null : gopts.palette;

  if (width <= 0 || height <= 0 || width > 65535 || height > 65535)
    throw "Width/Height invalid."

  function check_palette_and_num_colors(palette) {
    var num_colors = palette.length;
    if (num_colors < 2 || num_colors > 256 ||  num_colors & (num_colors-1))
      throw "Invalid code/color length, must be power of 2 and 2 .. 256.";
    return num_colors;
  }

  // - Header.
  buf[p++] = 0x47; buf[p++] = 0x49; buf[p++] = 0x46;  // GIF
  buf[p++] = 0x38; buf[p++] = 0x39; buf[p++] = 0x61;  // 89a

  // Handling of Global Color Table (palette) and background index.
  var gp_num_colors_pow2 = 0;
  var background = 0;
  if (global_palette !== null) {
    var gp_num_colors = check_palette_and_num_colors(global_palette);
    while (gp_num_colors >>= 1) ++gp_num_colors_pow2;
    gp_num_colors = 1 << gp_num_colors_pow2;
    --gp_num_colors_pow2;
    if (gopts.background !== undefined) {
      background = gopts.background;
      if (background >= gp_num_colors) throw "Background index out of range.";
      // The GIF spec states that a background index of 0 should be ignored, so
      // this is probably a mistake and you really want to set it to another
      // slot in the palette.  But actually in the end most browsers, etc end
      // up ignoring this almost completely (including for dispose background).
      if (background === 0)
        throw "Background index explicitly passed as 0.";
    }
  }

  // - Logical Screen Descriptor.
  // NOTE(deanm): w/h apparently ignored by implementations, but set anyway.
  buf[p++] = width & 0xff; buf[p++] = width >> 8 & 0xff;
  buf[p++] = height & 0xff; buf[p++] = height >> 8 & 0xff;
  // NOTE: Indicates 0-bpp original color resolution (unused?).
  buf[p++] = (global_palette !== null ? 0x80 : 0) |  // Global Color Table Flag.
             gp_num_colors_pow2;  // NOTE: No sort flag (unused?).
  buf[p++] = background;  // Background Color Index.
  buf[p++] = 0;  // Pixel aspect ratio (unused?).

  // - Global Color Table
  if (global_palette !== null) {
    for (var i = 0, il = global_palette.length; i < il; ++i) {
      var rgb = global_palette[i];
      buf[p++] = rgb >> 16 & 0xff;
      buf[p++] = rgb >> 8 & 0xff;
      buf[p++] = rgb & 0xff;
    }
  }

  if (loop_count !== null) {  // Netscape block for looping.
    if (loop_count < 0 || loop_count > 65535)
      throw "Loop count invalid."
    // Extension code, label, and length.
    buf[p++] = 0x21; buf[p++] = 0xff; buf[p++] = 0x0b;
    // NETSCAPE2.0
    buf[p++] = 0x4e; buf[p++] = 0x45; buf[p++] = 0x54; buf[p++] = 0x53;
    buf[p++] = 0x43; buf[p++] = 0x41; buf[p++] = 0x50; buf[p++] = 0x45;
    buf[p++] = 0x32; buf[p++] = 0x2e; buf[p++] = 0x30;
    // Sub-block
    buf[p++] = 0x03; buf[p++] = 0x01;
    buf[p++] = loop_count & 0xff; buf[p++] = loop_count >> 8 & 0xff;
    buf[p++] = 0x00;  // Terminator.
  }


  var ended = false;

  this.addFrame = function(x, y, w, h, indexed_pixels, opts) {
    if (ended === true) { --p; ended = false; }  // Un-end.

    opts = opts === undefined ? { } : opts;

    // TODO(deanm): Bounds check x, y.  Do they need to be within the virtual
    // canvas width/height, I imagine?
    if (x < 0 || y < 0 || x > 65535 || y > 65535)
      throw "x/y invalid."

    if (w <= 0 || h <= 0 || w > 65535 || h > 65535)
      throw "Width/Height invalid."

    if (indexed_pixels.length < w * h)
      throw "Not enough pixels for the frame size.";

    var using_local_palette = true;
    var palette = opts.palette;
    if (palette === undefined || palette === null) {
      using_local_palette = false;
      palette = global_palette;
    }

    if (palette === undefined || palette === null)
      throw "Must supply either a local or global palette.";

    var num_colors = check_palette_and_num_colors(palette);

    // Compute the min_code_size (power of 2), destroying num_colors.
    var min_code_size = 0;
    while (num_colors >>= 1) ++min_code_size;
    num_colors = 1 << min_code_size;  // Now we can easily get it back.

    var delay = opts.delay === undefined ? 0 : opts.delay;

    // From the spec:
    //     0 -   No disposal specified. The decoder is
    //           not required to take any action.
    //     1 -   Do not dispose. The graphic is to be left
    //           in place.
    //     2 -   Restore to background color. The area used by the
    //           graphic must be restored to the background color.
    //     3 -   Restore to previous. The decoder is required to
    //           restore the area overwritten by the graphic with
    //           what was there prior to rendering the graphic.
    //  4-7 -    To be defined.
    // NOTE(deanm): Dispose background doesn't really work, apparently most
    // browsers ignore the background palette index and clear to transparency.
    var disposal = opts.disposal === undefined ? 0 : opts.disposal;
    if (disposal < 0 || disposal > 3)  // 4-7 is reserved.
      throw "Disposal out of range.";

    var use_transparency = false;
    var transparent_index = 0;
    if (opts.transparent !== undefined && opts.transparent !== null) {
      use_transparency = true;
      transparent_index = opts.transparent;
      if (transparent_index < 0 || transparent_index >= num_colors)
        throw "Transparent color index.";
    }

    if (disposal !== 0 || use_transparency || delay !== 0) {
      // - Graphics Control Extension
      buf[p++] = 0x21; buf[p++] = 0xf9;  // Extension / Label.
      buf[p++] = 4;  // Byte size.

      buf[p++] = disposal << 2 | (use_transparency === true ? 1 : 0);
      buf[p++] = delay & 0xff; buf[p++] = delay >> 8 & 0xff;
      buf[p++] = transparent_index;  // Transparent color index.
      buf[p++] = 0;  // Block Terminator.
    }

    // - Image Descriptor
    buf[p++] = 0x2c;  // Image Seperator.
    buf[p++] = x & 0xff; buf[p++] = x >> 8 & 0xff;  // Left.
    buf[p++] = y & 0xff; buf[p++] = y >> 8 & 0xff;  // Top.
    buf[p++] = w & 0xff; buf[p++] = w >> 8 & 0xff;
    buf[p++] = h & 0xff; buf[p++] = h >> 8 & 0xff;
    // NOTE: No sort flag (unused?).
    // TODO(deanm): Support interlace.
    buf[p++] = using_local_palette === true ? (0x80 | (min_code_size-1)) : 0;

    // - Local Color Table
    if (using_local_palette === true) {
      for (var i = 0, il = palette.length; i < il; ++i) {
        var rgb = palette[i];
        buf[p++] = rgb >> 16 & 0xff;
        buf[p++] = rgb >> 8 & 0xff;
        buf[p++] = rgb & 0xff;
      }
    }

    p = GifWriterOutputLZWCodeStream(
            buf, p, min_code_size < 2 ? 2 : min_code_size, indexed_pixels);
  };

  this.end = function() {
    if (ended === false) {
      buf[p++] = 0x3b;  // Trailer.
      ended = true;
    }
    return p;
  };
}

// Main compression routine, palette indexes -> LZW code stream.
// |index_stream| must have at least one entry.
function GifWriterOutputLZWCodeStream(buf, p, min_code_size, index_stream) {
  buf[p++] = min_code_size;
  var cur_subblock = p++;  // Pointing at the length field.

  var clear_code = 1 << min_code_size;
  var code_mask = clear_code - 1;
  var eoi_code = clear_code + 1;
  var next_code = eoi_code + 1;

  var cur_code_size = min_code_size + 1;  // Number of bits per code.
  var cur_shift = 0;
  // We have at most 12-bit codes, so we should have to hold a max of 19
  // bits here (and then we would write out).
  var cur = 0;

  function emit_bytes_to_buffer(bit_block_size) {
    while (cur_shift >= bit_block_size) {
      buf[p++] = cur & 0xff;
      cur >>= 8; cur_shift -= 8;
      if (p === cur_subblock + 256) {  // Finished a subblock.
        buf[cur_subblock] = 255;
        cur_subblock = p++;
      }
    }
  }

  function emit_code(c) {
    cur |= c << cur_shift;
    cur_shift += cur_code_size;
    emit_bytes_to_buffer(8);
  }

  // I am not an expert on the topic, and I don't want to write a thesis.
  // However, it is good to outline here the basic algorithm and the few data
  // structures and optimizations here that make this implementation fast.
  // The basic idea behind LZW is to build a table of previously seen runs
  // addressed by a short id (herein called output code).  All data is
  // referenced by a code, which represents one or more values from the
  // original input stream.  All input bytes can be referenced as the same
  // value as an output code.  So if you didn't want any compression, you
  // could more or less just output the original bytes as codes (there are
  // some details to this, but it is the idea).  In order to achieve
  // compression, values greater then the input range (codes can be up to
  // 12-bit while input only 8-bit) represent a sequence of previously seen
  // inputs.  The decompressor is able to build the same mapping while
  // decoding, so there is always a shared common knowledge between the
  // encoding and decoder, which is also important for "timing" aspects like
  // how to handle variable bit width code encoding.
  //
  // One obvious but very important consequence of the table system is there
  // is always a unique id (at most 12-bits) to map the runs.  'A' might be
  // 4, then 'AA' might be 10, 'AAA' 11, 'AAAA' 12, etc.  This relationship
  // can be used for an effecient lookup strategy for the code mapping.  We
  // need to know if a run has been seen before, and be able to map that run
  // to the output code.  Since we start with known unique ids (input bytes),
  // and then from those build more unique ids (table entries), we can
  // continue this chain (almost like a linked list) to always have small
  // integer values that represent the current byte chains in the encoder.
  // This means instead of tracking the input bytes (AAAABCD) to know our
  // current state, we can track the table entry for AAAABC (it is guaranteed
  // to exist by the nature of the algorithm) and the next character D.
  // Therefor the tuple of (table_entry, byte) is guaranteed to also be
  // unique.  This allows us to create a simple lookup key for mapping input
  // sequences to codes (table indices) without having to store or search
  // any of the code sequences.  So if 'AAAA' has a table entry of 12, the
  // tuple of ('AAAA', K) for any input byte K will be unique, and can be our
  // key.  This leads to a integer value at most 20-bits, which can always
  // fit in an SMI value and be used as a fast sparse array / object key.

  // Output code for the current contents of the index buffer.
  var ib_code = index_stream[0] & code_mask;  // Load first input index.
  var code_table = { };  // Key'd on our 20-bit "tuple".

  emit_code(clear_code);  // Spec says first code should be a clear code.

  // First index already loaded, process the rest of the stream.
  for (var i = 1, il = index_stream.length; i < il; ++i) {
    var k = index_stream[i] & code_mask;
    var cur_key = ib_code << 8 | k;  // (prev, k) unique tuple.
    var cur_code = code_table[cur_key];  // buffer + k.

    // Check if we have to create a new code table entry.
    if (cur_code === undefined) {  // We don't have buffer + k.
      // Emit index buffer (without k).
      // This is an inline version of emit_code, because this is the core
      // writing routine of the compressor (and V8 cannot inline emit_code
      // because it is a closure here in a different context).  Additionally
      // we can call emit_byte_to_buffer less often, because we can have
      // 30-bits (from our 31-bit signed SMI), and we know our codes will only
      // be 12-bits, so can safely have 18-bits there without overflow.
      // emit_code(ib_code);
      cur |= ib_code << cur_shift;
      cur_shift += cur_code_size;
      while (cur_shift >= 8) {
        buf[p++] = cur & 0xff;
        cur >>= 8; cur_shift -= 8;
        if (p === cur_subblock + 256) {  // Finished a subblock.
          buf[cur_subblock] = 255;
          cur_subblock = p++;
        }
      }

      if (next_code === 4096) {  // Table full, need a clear.
        emit_code(clear_code);
        next_code = eoi_code + 1;
        cur_code_size = min_code_size + 1;
        code_table = { };
      } else {  // Table not full, insert a new entry.
        // Increase our variable bit code sizes if necessary.  This is a bit
        // tricky as it is based on "timing" between the encoding and
        // decoder.  From the encoders perspective this should happen after
        // we've already emitted the index buffer and are about to create the
        // first table entry that would overflow our current code bit size.
        if (next_code >= (1 << cur_code_size)) ++cur_code_size;
        code_table[cur_key] = next_code++;  // Insert into code table.
      }

      ib_code = k;  // Index buffer to single input k.
    } else {
      ib_code = cur_code;  // Index buffer to sequence in code table.
    }
  }

  emit_code(ib_code);  // There will still be something in the index buffer.
  emit_code(eoi_code);  // End Of Information.

  // Flush / finalize the sub-blocks stream to the buffer.
  emit_bytes_to_buffer(1);

  // Finish the sub-blocks, writing out any unfinished lengths and
  // terminating with a sub-block of length 0.  If we have already started
  // but not yet used a sub-block it can just become the terminator.
  if (cur_subblock + 1 === p) {  // Started but unused.
    buf[cur_subblock] = 0;
  } else {  // Started and used, write length and additional terminator block.
    buf[cur_subblock] = p - cur_subblock - 1;
    buf[p++] = 0;
  }
  return p;
}

function GifReader(buf) {
  var p = 0;

  // - Header (GIF87a or GIF89a).
  if (buf[p++] !== 0x47 ||            buf[p++] !== 0x49 || buf[p++] !== 0x46 ||
      buf[p++] !== 0x38 || (buf[p++]+1 & 0xfd) !== 0x38 || buf[p++] !== 0x61) {
    throw "Invalid GIF 87a/89a header.";
  }

  // - Logical Screen Descriptor.
  var width = buf[p++] | buf[p++] << 8;
  var height = buf[p++] | buf[p++] << 8;
  var pf0 = buf[p++];  // <Packed Fields>.
  var global_palette_flag = pf0 >> 7;
  var num_global_colors_pow2 = pf0 & 0x7;
  var num_global_colors = 1 << (num_global_colors_pow2 + 1);
  var background = buf[p++];
  buf[p++];  // Pixel aspect ratio (unused?).

  var global_palette_offset = null;

  if (global_palette_flag) {
    global_palette_offset = p;
    p += num_global_colors * 3;  // Seek past palette.
  }

  var no_eof = true;

  var frames = [ ];

  var delay = 0;
  var transparent_index = null;
  var disposal = 0;  // 0 - No disposal specified.
  var loop_count = null;

  this.width = width;
  this.height = height;

  while (no_eof && p < buf.length) {
    switch (buf[p++]) {
      case 0x21:  // Graphics Control Extension Block
        switch (buf[p++]) {
          case 0xff:  // Application specific block
            // Try if it's a Netscape block (with animation loop counter).
            if (buf[p   ] !== 0x0b ||  // 21 FF already read, check block size.
                // NETSCAPE2.0
                buf[p+1 ] == 0x4e && buf[p+2 ] == 0x45 && buf[p+3 ] == 0x54 &&
                buf[p+4 ] == 0x53 && buf[p+5 ] == 0x43 && buf[p+6 ] == 0x41 &&
                buf[p+7 ] == 0x50 && buf[p+8 ] == 0x45 && buf[p+9 ] == 0x32 &&
                buf[p+10] == 0x2e && buf[p+11] == 0x30 &&
                // Sub-block
                buf[p+12] == 0x03 && buf[p+13] == 0x01 && buf[p+16] == 0) {
              p += 14;
              loop_count = buf[p++] | buf[p++] << 8;
              p++;  // Skip terminator.
            } else {  // We don't know what it is, just try to get past it.
              p += 12;
              while (true) {  // Seek through subblocks.
                var block_size = buf[p++];
                if (block_size === 0) break;
                p += block_size;
              }
            }
            break;

          case 0xf9:  // Graphics Control Extension
            if (buf[p++] !== 0x4 || buf[p+4] !== 0)
              throw "Invalid graphics extension block.";
            var pf1 = buf[p++];
            delay = buf[p++] | buf[p++] << 8;
            transparent_index = buf[p++];
            if ((pf1 & 1) === 0) transparent_index = null;
            disposal = pf1 >> 2 & 0x7;
            p++;  // Skip terminator.
            break;

          case 0xfe:  // Comment Extension.
            while (true) {  // Seek through subblocks.
              var block_size = buf[p++];
              if (block_size === 0) break;
              // console.log(buf.slice(p, p+block_size).toString('ascii'));
              p += block_size;
            }
            break;

          default:
            throw "Unknown graphic control label: 0x" + buf[p-1].toString(16);
        }
        break;

      case 0x2c:  // Image Descriptor.
        var x = buf[p++] | buf[p++] << 8;
        var y = buf[p++] | buf[p++] << 8;
        var w = buf[p++] | buf[p++] << 8;
        var h = buf[p++] | buf[p++] << 8;
        var pf2 = buf[p++];
        var local_palette_flag = pf2 >> 7;
        var interlace_flag = pf2 >> 6 & 1;
        var num_local_colors_pow2 = pf2 & 0x7;
        var num_local_colors = 1 << (num_local_colors_pow2 + 1);
        var palette_offset = global_palette_offset;
        var has_local_palette = false;
        if (local_palette_flag) {
          var has_local_palette = true;
          palette_offset = p;  // Override with local palette.
          p += num_local_colors * 3;  // Seek past palette.
        }

        var data_offset = p;

        p++;  // codesize
        while (true) {
          var block_size = buf[p++];
          if (block_size === 0) break;
          p += block_size;
        }

        frames.push({x: x, y: y, width: w, height: h,
                     has_local_palette: has_local_palette,
                     palette_offset: palette_offset,
                     data_offset: data_offset,
                     data_length: p - data_offset,
                     transparent_index: transparent_index,
                     interlaced: !!interlace_flag,
                     delay: delay,
                     disposal: disposal});
        break;

      case 0x3b:  // Trailer Marker (end of file).
        no_eof = false;
        break;

      default:
        throw "Unknown gif block: 0x" + buf[p-1].toString(16);
        break;
    }
  }

  this.numFrames = function() {
    return frames.length;
  };

  this.loopCount = function() {
    return loop_count;
  };

  this.frameInfo = function(frame_num) {
    if (frame_num < 0 || frame_num >= frames.length)
      throw "Frame index out of range.";
    return frames[frame_num];
  }

  this.decodeAndBlitFrameBGRA = function(frame_num, pixels) {
    var frame = this.frameInfo(frame_num);
    var num_pixels = frame.width * frame.height;
    var index_stream = new Uint8Array(num_pixels);  // At most 8-bit indices.
    GifReaderLZWOutputIndexStream(
        buf, frame.data_offset, index_stream, num_pixels);
    var palette_offset = frame.palette_offset;

    // NOTE(deanm): It seems to be much faster to compare index to 256 than
    // to === null.  Not sure why, but CompareStub_EQ_STRICT shows up high in
    // the profile, not sure if it's related to using a Uint8Array.
    var trans = frame.transparent_index;
    if (trans === null) trans = 256;

    // We are possibly just blitting to a portion of the entire frame.
    // That is a subrect within the framerect, so the additional pixels
    // must be skipped over after we finished a scanline.
    var framewidth  = frame.width;
    var framestride = width - framewidth;
    var xleft       = framewidth;  // Number of subrect pixels left in scanline.

    // Output indicies of the top left and bottom right corners of the subrect.
    var opbeg = ((frame.y * width) + frame.x) * 4;
    var opend = ((frame.y + frame.height) * width + frame.x) * 4;
    var op    = opbeg;

    var scanstride = framestride * 4;

    // Use scanstride to skip past the rows when interlacing.  This is skipping
    // 7 rows for the first two passes, then 3 then 1.
    if (frame.interlaced === true) {
      scanstride += width * 4 * 7;  // Pass 1.
    }

    var interlaceskip = 8;  // Tracking the row interval in the current pass.

    for (var i = 0, il = index_stream.length; i < il; ++i) {
      var index = index_stream[i];

      if (xleft === 0) {  // Beginning of new scan line
        op += scanstride;
        xleft = framewidth;
        if (op >= opend) { // Catch the wrap to switch passes when interlacing.
          scanstride = framestride * 4 + width * 4 * (interlaceskip-1);
          // interlaceskip / 2 * 4 is interlaceskip << 1.
          op = opbeg + (framewidth + framestride) * (interlaceskip << 1);
          interlaceskip >>= 1;
        }
      }

      if (index === trans) {
        op += 4;
      } else {
        var r = buf[palette_offset + index * 3];
        var g = buf[palette_offset + index * 3 + 1];
        var b = buf[palette_offset + index * 3 + 2];
        pixels[op++] = b;
        pixels[op++] = g;
        pixels[op++] = r;
        pixels[op++] = 255;
      }
      --xleft;
    }
  };

  // I will go to copy and paste hell one day...
  this.decodeAndBlitFrameRGBA = function(frame_num, pixels) {
    var frame = this.frameInfo(frame_num);
    var num_pixels = frame.width * frame.height;
    var index_stream = new Uint8Array(num_pixels);  // At most 8-bit indices.
    GifReaderLZWOutputIndexStream(
        buf, frame.data_offset, index_stream, num_pixels);
    var palette_offset = frame.palette_offset;

    // NOTE(deanm): It seems to be much faster to compare index to 256 than
    // to === null.  Not sure why, but CompareStub_EQ_STRICT shows up high in
    // the profile, not sure if it's related to using a Uint8Array.
    var trans = frame.transparent_index;
    if (trans === null) trans = 256;

    // We are possibly just blitting to a portion of the entire frame.
    // That is a subrect within the framerect, so the additional pixels
    // must be skipped over after we finished a scanline.
    var framewidth  = frame.width;
    var framestride = width - framewidth;
    var xleft       = framewidth;  // Number of subrect pixels left in scanline.

    // Output indicies of the top left and bottom right corners of the subrect.
    var opbeg = ((frame.y * width) + frame.x) * 4;
    var opend = ((frame.y + frame.height) * width + frame.x) * 4;
    var op    = opbeg;

    var scanstride = framestride * 4;

    // Use scanstride to skip past the rows when interlacing.  This is skipping
    // 7 rows for the first two passes, then 3 then 1.
    if (frame.interlaced === true) {
      scanstride += width * 4 * 7;  // Pass 1.
    }

    var interlaceskip = 8;  // Tracking the row interval in the current pass.

    for (var i = 0, il = index_stream.length; i < il; ++i) {
      var index = index_stream[i];

      if (xleft === 0) {  // Beginning of new scan line
        op += scanstride;
        xleft = framewidth;
        if (op >= opend) { // Catch the wrap to switch passes when interlacing.
          scanstride = framestride * 4 + width * 4 * (interlaceskip-1);
          // interlaceskip / 2 * 4 is interlaceskip << 1.
          op = opbeg + (framewidth + framestride) * (interlaceskip << 1);
          interlaceskip >>= 1;
        }
      }

      if (index === trans) {
        op += 4;
      } else {
        var r = buf[palette_offset + index * 3];
        var g = buf[palette_offset + index * 3 + 1];
        var b = buf[palette_offset + index * 3 + 2];
        pixels[op++] = r;
        pixels[op++] = g;
        pixels[op++] = b;
        pixels[op++] = 255;
      }
      --xleft;
    }
  };
}

function GifReaderLZWOutputIndexStream(code_stream, p, output, output_length) {
  var min_code_size = code_stream[p++];

  var clear_code = 1 << min_code_size;
  var eoi_code = clear_code + 1;
  var next_code = eoi_code + 1;

  var cur_code_size = min_code_size + 1;  // Number of bits per code.
  // NOTE: This shares the same name as the encoder, but has a different
  // meaning here.  Here this masks each code coming from the code stream.
  var code_mask = (1 << cur_code_size) - 1;
  var cur_shift = 0;
  var cur = 0;

  var op = 0;  // Output pointer.
  
  var subblock_size = code_stream[p++];

  // TODO(deanm): Would using a TypedArray be any faster?  At least it would
  // solve the fast mode / backing store uncertainty.
  // var code_table = Array(4096);
  var code_table = new Int32Array(4096);  // Can be signed, we only use 20 bits.

  var prev_code = null;  // Track code-1.

  while (true) {
    // Read up to two bytes, making sure we always 12-bits for max sized code.
    while (cur_shift < 16) {
      if (subblock_size === 0) break;  // No more data to be read.

      cur |= code_stream[p++] << cur_shift;
      cur_shift += 8;

      if (subblock_size === 1) {  // Never let it get to 0 to hold logic above.
        subblock_size = code_stream[p++];  // Next subblock.
      } else {
        --subblock_size;
      }
    }

    // TODO(deanm): We should never really get here, we should have received
    // and EOI.
    if (cur_shift < cur_code_size)
      break;

    var code = cur & code_mask;
    cur >>= cur_code_size;
    cur_shift -= cur_code_size;

    // TODO(deanm): Maybe should check that the first code was a clear code,
    // at least this is what you're supposed to do.  But actually our encoder
    // now doesn't emit a clear code first anyway.
    if (code === clear_code) {
      // We don't actually have to clear the table.  This could be a good idea
      // for greater error checking, but we don't really do any anyway.  We
      // will just track it with next_code and overwrite old entries.

      next_code = eoi_code + 1;
      cur_code_size = min_code_size + 1;
      code_mask = (1 << cur_code_size) - 1;

      // Don't update prev_code ?
      prev_code = null;
      continue;
    } else if (code === eoi_code) {
      break;
    }

    // We have a similar situation as the decoder, where we want to store
    // variable length entries (code table entries), but we want to do in a
    // faster manner than an array of arrays.  The code below stores sort of a
    // linked list within the code table, and then "chases" through it to
    // construct the dictionary entries.  When a new entry is created, just the
    // last byte is stored, and the rest (prefix) of the entry is only
    // referenced by its table entry.  Then the code chases through the
    // prefixes until it reaches a single byte code.  We have to chase twice,
    // first to compute the length, and then to actually copy the data to the
    // output (backwards, since we know the length).  The alternative would be
    // storing something in an intermediate stack, but that doesn't make any
    // more sense.  I implemented an approach where it also stored the length
    // in the code table, although it's a bit tricky because you run out of
    // bits (12 + 12 + 8), but I didn't measure much improvements (the table
    // entries are generally not the long).  Even when I created benchmarks for
    // very long table entries the complexity did not seem worth it.
    // The code table stores the prefix entry in 12 bits and then the suffix
    // byte in 8 bits, so each entry is 20 bits.

    var chase_code = code < next_code ? code : prev_code;

    // Chase what we will output, either {CODE} or {CODE-1}.
    var chase_length = 0;
    var chase = chase_code;
    while (chase > clear_code) {
      chase = code_table[chase] >> 8;
      ++chase_length;
    }

    var k = chase;
    
    var op_end = op + chase_length + (chase_code !== code ? 1 : 0);
    if (op_end > output_length) {
      console.log("Warning, gif stream longer than expected.");
      return;
    }

    // Already have the first byte from the chase, might as well write it fast.
    output[op++] = k;

    op += chase_length;
    var b = op;  // Track pointer, writing backwards.

    if (chase_code !== code)  // The case of emitting {CODE-1} + k.
      output[op++] = k;

    chase = chase_code;
    while (chase_length--) {
      chase = code_table[chase];
      output[--b] = chase & 0xff;  // Write backwards.
      chase >>= 8;  // Pull down to the prefix code.
    }

    if (prev_code !== null && next_code < 4096) {
      code_table[next_code++] = prev_code << 8 | k;
      // TODO(deanm): Figure out this clearing vs code growth logic better.  I
      // have an feeling that it should just happen somewhere else, for now it
      // is awkward between when we grow past the max and then hit a clear code.
      // For now just check if we hit the max 12-bits (then a clear code should
      // follow, also of course encoded in 12-bits).
      if (next_code >= code_mask+1 && cur_code_size < 12) {
        ++cur_code_size;
        code_mask = code_mask << 1 | 1;
      }
    }

    prev_code = code;
  }

  if (op !== output_length) {
    console.log("Warning, gif stream shorter than expected.");
  }

  return output;
}

try { exports.GifWriter = GifWriter; exports.GifReader = GifReader } catch(e) { }  // CommonJS.

},{}],4:[function(require,module,exports){
var GifReader;

GifReader = require('omggif').GifReader;

Phaser.Loader.prototype.spritegif = function(key, url) {
  console.log('GIF: ', key, url);
  return this.binary(key, url, ((function(_this) {
    return function(key, data) {
      var canvas, context, i, image, index, reader, ref, ref1, square, start, value, x, y;
      start = Date.now();
      reader = new GifReader(new Uint8Array(data));
      square = Math.ceil(Math.sqrt(reader.numFrames()));
      console.log(square * reader.width, square * reader.height);
      canvas = Phaser.Canvas.create(square * reader.width, square * reader.height);
      context = canvas.getContext('2d');
      for (index = i = 0, ref = reader.numFrames(); 0 <= ref ? i < ref : i > ref; index = 0 <= ref ? ++i : --i) {
        x = (index % square) * reader.width;
        y = (Math.floor(index / square)) * reader.height;
        image = new ImageData(reader.width, reader.height);
        ref1 = reader.frameInfo(index);
        for (key in ref1) {
          value = ref1[key];
          image[key] = value;
        }
        if (image.delay) {
          image.delay = image.delay * 10;
        }
        reader.decodeAndBlitFrameRGBA(index, image.data);
        context.putImageData(image, x, y);
      }
      return Phaser.Canvas.addToDOM(canvas, document.body);
    };
  })(this)), this.game);
};


},{"omggif":3}],5:[function(require,module,exports){
var Luna, lune,
  extend = function(child, parent) { for (var key in parent) { if (hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; },
  hasProp = {}.hasOwnProperty;

lune = require('lune');

console.log(lune.phase());

require('./gif.coffee');

Luna = (function(superClass) {
  extend(Luna, superClass);

  function Luna() {
    if (!(this instanceof Luna)) {
      return new Luna;
    }
    this.game = Luna.__super__.constructor.call(this, window.innerWidth, window.innerHeight, Phaser.AUTO, '', this);
  }

  Luna.prototype.preload = function() {
    this.game.load.crossOrigin = '*';
    return this.game.load.spritegif('m8hYEjb', 'http://i.imgur.com/m8hYEjb.gif');
  };

  Luna.prototype.create = function() {};

  return Luna;

})(Phaser.Game);

window.game = new Luna;


},{"./gif.coffee":4,"lune":1}]},{},[5])
//# sourceMappingURL=data:application/json;charset:utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIm5vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCJub2RlX21vZHVsZXMvbHVuZS9pbmRleC5qcyIsIm5vZGVfbW9kdWxlcy9sdW5lL2xpYi9sdW5lLmpzIiwibm9kZV9tb2R1bGVzL29tZ2dpZi9vbWdnaWYuanMiLCJDOlxcZGV2XFxub2RlX21vZHVsZXNcXHByb2RcXG5vZGVfbW9kdWxlc1xcbHVuYVxcc3JjXFxnaWYuY29mZmVlIiwiQzpcXGRldlxcbm9kZV9tb2R1bGVzXFxwcm9kXFxub2RlX21vZHVsZXNcXGx1bmFcXHNyY1xcbHVuYS5jb2ZmZWUiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBOztBQ0RBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3JjQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDOXdCQSxJQUFBOztBQUFDLFlBQWEsT0FBQSxDQUFRLFFBQVIsRUFBYjs7QUFFRCxNQUFNLENBQUMsTUFBTSxDQUFDLFNBQVMsQ0FBQyxTQUF4QixHQUFvQyxTQUFDLEdBQUQsRUFBTSxHQUFOO0VBQ2xDLE9BQU8sQ0FBQyxHQUFSLENBQVksT0FBWixFQUFxQixHQUFyQixFQUEwQixHQUExQjtTQUNBLElBQUMsQ0FBQSxNQUFELENBQVEsR0FBUixFQUFhLEdBQWIsRUFBa0IsQ0FDaEIsQ0FBQSxTQUFBLEtBQUE7V0FBQSxTQUFDLEdBQUQsRUFBTSxJQUFOO0FBQ0UsVUFBQTtNQUFBLEtBQUEsR0FBUSxJQUFJLENBQUMsR0FBTCxDQUFBO01BQ1IsTUFBQSxHQUFhLElBQUEsU0FBQSxDQUFjLElBQUEsVUFBQSxDQUFXLElBQVgsQ0FBZDtNQUNiLE1BQUEsR0FBUyxJQUFJLENBQUMsSUFBTCxDQUFVLElBQUksQ0FBQyxJQUFMLENBQVUsTUFBTSxDQUFDLFNBQVAsQ0FBQSxDQUFWLENBQVY7TUFDVCxPQUFPLENBQUMsR0FBUixDQUFZLE1BQUEsR0FBUyxNQUFNLENBQUMsS0FBNUIsRUFBbUMsTUFBQSxHQUFTLE1BQU0sQ0FBQyxNQUFuRDtNQUNBLE1BQUEsR0FBUyxNQUFNLENBQUMsTUFBTSxDQUFDLE1BQWQsQ0FBcUIsTUFBQSxHQUFTLE1BQU0sQ0FBQyxLQUFyQyxFQUE0QyxNQUFBLEdBQVMsTUFBTSxDQUFDLE1BQTVEO01BQ1QsT0FBQSxHQUFVLE1BQU0sQ0FBQyxVQUFQLENBQWtCLElBQWxCO0FBQ1YsV0FBYSxtR0FBYjtRQUNFLENBQUEsR0FBSSxDQUFDLEtBQUEsR0FBUSxNQUFULENBQUEsR0FBbUIsTUFBTSxDQUFDO1FBQzlCLENBQUEsR0FBSSxDQUFDLElBQUksQ0FBQyxLQUFMLENBQVcsS0FBQSxHQUFRLE1BQW5CLENBQUQsQ0FBQSxHQUE4QixNQUFNLENBQUM7UUFDekMsS0FBQSxHQUFZLElBQUEsU0FBQSxDQUFVLE1BQU0sQ0FBQyxLQUFqQixFQUF3QixNQUFNLENBQUMsTUFBL0I7QUFDWjtBQUFBLGFBQUEsV0FBQTs7VUFBQSxLQUFNLENBQUEsR0FBQSxDQUFOLEdBQWE7QUFBYjtRQUNBLElBQWtDLEtBQUssQ0FBQyxLQUF4QztVQUFBLEtBQUssQ0FBQyxLQUFOLEdBQWMsS0FBSyxDQUFDLEtBQU4sR0FBYyxHQUE1Qjs7UUFDQSxNQUFNLENBQUMsc0JBQVAsQ0FBOEIsS0FBOUIsRUFBcUMsS0FBSyxDQUFDLElBQTNDO1FBQ0EsT0FBTyxDQUFDLFlBQVIsQ0FBcUIsS0FBckIsRUFBNEIsQ0FBNUIsRUFBK0IsQ0FBL0I7QUFQRjthQVFBLE1BQU0sQ0FBQyxNQUFNLENBQUMsUUFBZCxDQUF1QixNQUF2QixFQUErQixRQUFRLENBQUMsSUFBeEM7SUFmRjtFQUFBLENBQUEsQ0FBQSxDQUFBLElBQUEsQ0FEZ0IsQ0FBbEIsRUFpQkcsSUFBQyxDQUFBLElBakJKO0FBRmtDOzs7O0FDRnBDLElBQUEsVUFBQTtFQUFBOzs7QUFBQSxJQUFBLEdBQU8sT0FBQSxDQUFRLE1BQVI7O0FBQ1AsT0FBTyxDQUFDLEdBQVIsQ0FBWSxJQUFJLENBQUMsS0FBTCxDQUFBLENBQVo7O0FBRUEsT0FBQSxDQUFRLGNBQVI7O0FBRU07OztFQUNTLGNBQUE7SUFDWCxJQUFBLENBQUEsQ0FBTyxJQUFBLFlBQWEsSUFBcEIsQ0FBQTtBQUE4QixhQUFPLElBQUksS0FBekM7O0lBQ0EsSUFBQyxDQUFBLElBQUQsR0FBUSxzQ0FBTSxNQUFNLENBQUMsVUFBYixFQUF5QixNQUFNLENBQUMsV0FBaEMsRUFBNkMsTUFBTSxDQUFDLElBQXBELEVBQTBELEVBQTFELEVBQThELElBQTlEO0VBRkc7O2lCQUdiLE9BQUEsR0FBUyxTQUFBO0lBQ1AsSUFBQyxDQUFBLElBQUksQ0FBQyxJQUFJLENBQUMsV0FBWCxHQUF5QjtXQUN6QixJQUFDLENBQUEsSUFBSSxDQUFDLElBQUksQ0FBQyxTQUFYLENBQXFCLFNBQXJCLEVBQWdDLGdDQUFoQztFQUZPOztpQkFJVCxNQUFBLEdBQVEsU0FBQSxHQUFBOzs7O0dBUlMsTUFBTSxDQUFDOztBQVUxQixNQUFNLENBQUMsSUFBUCxHQUFjLElBQUkiLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbiBlKHQsbixyKXtmdW5jdGlvbiBzKG8sdSl7aWYoIW5bb10pe2lmKCF0W29dKXt2YXIgYT10eXBlb2YgcmVxdWlyZT09XCJmdW5jdGlvblwiJiZyZXF1aXJlO2lmKCF1JiZhKXJldHVybiBhKG8sITApO2lmKGkpcmV0dXJuIGkobywhMCk7dmFyIGY9bmV3IEVycm9yKFwiQ2Fubm90IGZpbmQgbW9kdWxlICdcIitvK1wiJ1wiKTt0aHJvdyBmLmNvZGU9XCJNT0RVTEVfTk9UX0ZPVU5EXCIsZn12YXIgbD1uW29dPXtleHBvcnRzOnt9fTt0W29dWzBdLmNhbGwobC5leHBvcnRzLGZ1bmN0aW9uKGUpe3ZhciBuPXRbb11bMV1bZV07cmV0dXJuIHMobj9uOmUpfSxsLGwuZXhwb3J0cyxlLHQsbixyKX1yZXR1cm4gbltvXS5leHBvcnRzfXZhciBpPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7Zm9yKHZhciBvPTA7bzxyLmxlbmd0aDtvKyspcyhyW29dKTtyZXR1cm4gc30pIiwibW9kdWxlLmV4cG9ydHMgPSByZXF1aXJlKCcuL2xpYi9sdW5lLmpzJyk7XG4iLCIvKipcbiAqIFRoaXMgbGlicmFyeSBjYWxjdWxhdGVzIHRoZSBjdXJyZW50IHBoYXNlIG9mIHRoZSBtb29uXG4gKiBhcyB3ZWxsIGFzIGZpbmRzIHRoZSBkYXRlcyBvZiB0aGUgcmVjZW50IG1vb24gcGhhc2VzLlxuICpcbiAqIFBvcnRlZCBmcm9tIHB5dGhvbiB2ZXJzaW9uIGZvdW5kIGhlcmU6XG4gKiBodHRwczovL2JhemFhci5sYXVuY2hwYWQubmV0L35rZXR1cm4vcHktbW9vbi1waGFzZS90cnVuay9hbm5vdGF0ZS9oZWFkOi9tb29uLnB5XG4gKlxuICogQXV0aG9yOiBSeWFuIFNleXMgKGh0dHBzOi8vZ2l0aHViLmNvbS9yeWFuc2V5cylcbiAqL1xuXG4oZnVuY3Rpb24oKSB7XG4gIC8vIFBoYXNlcyBvZiB0aGUgbW9vbiAmIHByZWNpc2lvblxuXG4gIHZhciBQUkVDSVNJT04gPSAwLjA1O1xuICB2YXIgTkVXID0gMCAvIDQuMDtcbiAgdmFyIEZJUlNUID0gMSAvIDQuMDtcbiAgdmFyIEZVTEwgPSAyIC8gNC4wO1xuICB2YXIgTEFTVCA9IDMgLyA0LjA7XG4gIHZhciBORVhUTkVXID0gNCAvIDQuMDtcblxuICAvKipcbiAgICogR2V0cyB0aGUgSnVsaWFuIHZhbHVlIGZyb20gYSBkYXRlIG9iamVjdC5cbiAgICogU291cmNlOiBodHRwOi8vamF2YXNjcmlwdC5hYm91dC5jb20vbGlicmFyeS9ibGp1bGRheS5odG1cbiAgICogQHJldHVybiB7TnVtYmVyfSBKdWxpYW4gbnVtYmVyIHJlcHJlc2VudGF0aW9uIG9mIHRoZSBkYXRlLlxuICAgKi9cbiAgRGF0ZS5wcm90b3R5cGUuZ2V0SnVsaWFuID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuICh0aGlzLnZhbHVlT2YoKSAvIDg2NDAwMDAwKSAtICh0aGlzLmdldFRpbWV6b25lT2Zmc2V0KCkgLyAxNDQwKSArIDI0NDA1ODcuNTtcbiAgfTtcblxuICAvKipcbiAgICogQ29udmVydHMgYSBOdW1iZXIgaW4gSnVsaWFuIGRhdGUgZm9ybSB0byBhIERhdGUgb2JqZWN0LlxuICAgKiBTb3VyY2U6IGh0dHA6Ly9ibG9nLmJhaHJlbmJ1cmdzLmNvbS8yMDExLzAxL2phdmFzY3JpcHQtanVsaWFuLWRheS1jb252ZXJzaW9ucy5odG1sXG4gICAqL1xuICBOdW1iZXIucHJvdG90eXBlLkp1bGlhbjJEYXRlID0gZnVuY3Rpb24oaW5VVEMpIHtcbiAgICB2YXIgWCA9IHBhcnNlRmxvYXQodGhpcykrMC41O1xuICAgIHZhciBaID0gTWF0aC5mbG9vcihYKTsgLy9HZXQgZGF5IHdpdGhvdXQgdGltZVxuICAgIHZhciBGID0gWCAtIFo7IC8vR2V0IHRpbWVcbiAgICB2YXIgWSA9IE1hdGguZmxvb3IoKFotMTg2NzIxNi4yNSkvMzY1MjQuMjUpO1xuICAgIHZhciBBID0gWisxK1ktTWF0aC5mbG9vcihZLzQpO1xuICAgIHZhciBCID0gQSsxNTI0O1xuICAgIHZhciBDID0gTWF0aC5mbG9vcigoQi0xMjIuMSkvMzY1LjI1KTtcbiAgICB2YXIgRCA9IE1hdGguZmxvb3IoMzY1LjI1KkMpO1xuICAgIHZhciBHID0gTWF0aC5mbG9vcigoQi1EKS8zMC42MDAxKTtcbiAgICAvL211c3QgZ2V0IG51bWJlciBsZXNzIHRoYW4gb3IgZXF1YWwgdG8gMTIpXG4gICAgdmFyIG1vbnRoID0gKEc8MTMuNSkgPyAoRy0xKSA6IChHLTEzKTtcbiAgICAvL2lmIE1vbnRoIGlzIEphbnVhcnkgb3IgRmVicnVhcnksIG9yIHRoZSByZXN0IG9mIHllYXJcbiAgICB2YXIgeWVhciA9IChtb250aDwyLjUpID8gKEMtNDcxNSkgOiAoQy00NzE2KTtcbiAgICBtb250aCAtPSAxOyAvL0hhbmRsZSBKYXZhU2NyaXB0IG1vbnRoIGZvcm1hdFxuICAgIHZhciBVVCA9IEItRC1NYXRoLmZsb29yKDMwLjYwMDEqRykrRjtcbiAgICB2YXIgZGF5ID0gTWF0aC5mbG9vcihVVCk7XG4gICAgLy9EZXRlcm1pbmUgdGltZVxuICAgIFVUIC09IE1hdGguZmxvb3IoVVQpO1xuICAgIFVUICo9IDI0O1xuICAgIHZhciBob3VyID0gTWF0aC5mbG9vcihVVCk7XG4gICAgVVQgLT0gTWF0aC5mbG9vcihVVCk7XG4gICAgVVQgKj0gNjA7XG4gICAgdmFyIG1pbnV0ZSA9IE1hdGguZmxvb3IoVVQpO1xuICAgIFVUIC09IE1hdGguZmxvb3IoVVQpO1xuICAgIFVUICo9IDYwO1xuICAgIHZhciBzZWNvbmQgPSBNYXRoLnJvdW5kKFVUKTtcblxuICAgIGlmIChpblVUQykge1xuICAgICAgcmV0dXJuIG5ldyBEYXRlKERhdGUuVVRDKHllYXIsIG1vbnRoLCBkYXksIGhvdXIsIG1pbnV0ZSwgc2Vjb25kKSk7XG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBuZXcgRGF0ZSh5ZWFyLCBtb250aCwgZGF5LCBob3VyLCBtaW51dGUsIHNlY29uZCk7XG4gICAgfVxuICB9O1xuXG4gIC8qKlxuICAgKiBBc3Ryb25vbWljYWwgQ29uc3RhbnRzXG4gICAqIEB0eXBlIHtPYmplY3R9XG4gICAqL1xuICBjb25zdCBjID0ge1xuICAgIC8vIEpETiBzdGFuZHMgZm9yIEp1bGlhbiBEYXkgTnVtYmVyXG4gICAgLy8gQW5nbGVzIGhlcmUgYXJlIGluIGRlZ3JlZXNcblxuICAgIC8vIDE5ODAgSmFudWFyeSAwLjAgaW4gSkROXG4gICAgLy8gWFhYOiBEYXRlVGltZSgxOTgwKS5qZG4geWllbGRzIDI0NDQyMzkuNSAtLSB3aGljaCBvbmUgaXMgcmlnaHQ/XG4gICAgZXBvY2g6IDI0NDQyMzguNSxcblxuICAgIC8vIEVjbGlwdGljIGxvbmdpdHVkZSBvZiB0aGUgU3VuIGF0IGVwb2NoIDE5ODAuMFxuICAgIGVjbGlwdGljX2xvbmdpdHVkZV9lcG9jaDogMjc4LjgzMzU0MCxcblxuICAgIC8vIEVjbGlwdGljIGxvbmdpdHVkZSBvZiB0aGUgU3VuIGF0IHBlcmlnZWVcbiAgICBlY2xpcHRpY19sb25naXR1ZGVfcGVyaWdlZTogMjgyLjU5NjQwMyxcblxuICAgIC8vIEVjY2VudHJpY2l0eSBvZiBFYXJ0aCdzIG9yYml0XG4gICAgZWNjZW50cmljaXR5OiAwLjAxNjcxOCxcblxuICAgIC8vIFNlbWktbWFqb3IgYXhpcyBvZiBFYXJ0aCdzIG9yYml0LCBpbiBraWxvbWV0ZXJzXG4gICAgc3VuX3NtYXhpczogMS40OTU4NWU4LFxuXG4gICAgLy8gU3VuJ3MgYW5ndWxhciBzaXplLCBpbiBkZWdyZWVzLCBhdCBzZW1pLW1ham9yIGF4aXMgZGlzdGFuY2VcbiAgICBzdW5fYW5ndWxhcl9zaXplX3NtYXhpczogMC41MzMxMjgsXG5cbiAgICAvLyBFbGVtZW50cyBvZiB0aGUgTW9vbidzIG9yYml0LCBlcG9jaCAxOTgwLjBcblxuICAgIC8vIE1vb24ncyBtZWFuIGxvbmdpdHVkZSBhdCB0aGUgZXBvY2hcbiAgICBtb29uX21lYW5fbG9uZ2l0dWRlX2Vwb2NoOiA2NC45NzU0NjQsXG4gICAgLy8gTWVhbiBsb25naXR1ZGUgb2YgdGhlIHBlcmlnZWUgYXQgdGhlIGVwb2NoXG4gICAgbW9vbl9tZWFuX3BlcmlnZWVfZXBvY2g6IDM0OS4zODMwNjMsXG5cbiAgICAvLyBNZWFuIGxvbmdpdHVkZSBvZiB0aGUgbm9kZSBhdCB0aGUgZXBvY2hcbiAgICBub2RlX21lYW5fbG9uZ2l0dWRlX2Vwb2NoOiAxNTEuOTUwNDI5LFxuXG4gICAgLy8gSW5jbGluYXRpb24gb2YgdGhlIE1vb24ncyBvcmJpdFxuICAgIG1vb25faW5jbGluYXRpb246IDUuMTQ1Mzk2LFxuXG4gICAgLy8gRWNjZW50cmljaXR5IG9mIHRoZSBNb29uJ3Mgb3JiaXRcbiAgICBtb29uX2VjY2VudHJpY2l0eTogMC4wNTQ5MDAsXG5cbiAgICAvLyBNb29uJ3MgYW5ndWxhciBzaXplIGF0IGRpc3RhbmNlIGEgZnJvbSBFYXJ0aFxuICAgIG1vb25fYW5ndWxhcl9zaXplOiAwLjUxODEsXG5cbiAgICAvLyBTZW1pLW1vam9yIGF4aXMgb2YgdGhlIE1vb24ncyBvcmJpdCwgaW4ga2lsb21ldGVyc1xuICAgIG1vb25fc21heGlzOiAzODQ0MDEuMCxcbiAgICAvLyBQYXJhbGxheCBhdCBhIGRpc3RhbmNlIGEgZnJvbSBFYXJ0aFxuICAgIG1vb25fcGFyYWxsYXg6IDAuOTUwNyxcblxuICAgIC8vIFN5bm9kaWMgbW9udGggKG5ldyBNb29uIHRvIG5ldyBNb29uKSwgaW4gZGF5c1xuICAgIHN5bm9kaWNfbW9udGg6IDI5LjUzMDU4ODY4LFxuXG4gICAgLy8gQmFzZSBkYXRlIGZvciBFLiBXLiBCcm93bidzIG51bWJlcmVkIHNlcmllcyBvZiBsdW5hdGlvbnMgKDE5MjMgSmFudWFyeSAxNilcbiAgICBsdW5hdGlvbnNfYmFzZTogMjQyMzQzNi4wLFxuXG4gICAgLy8gI1Byb3BlcnRpZXMgb2YgdGhlIEVhcnRoXG4gICAgZWFydGhfcmFkaXVzOiA2Mzc4LjE2XG4gIH07XG5cbiAgZnVuY3Rpb24gZml4YW5nbGUoYSkge1xuICAgIHJldHVybiBhIC0gMzYwLjAgKiBNYXRoLmZsb29yKGEvMzYwLjApO1xuICB9XG5cbiAgLyoqXG4gICAqIENvbnZlcnQgZGVncmVlcyB0byByYWRpYW5zXG4gICAqIEBwYXJhbSAge051bWJlcn0gZCBBbmdsZSBpbiBkZWdyZWVzXG4gICAqIEByZXR1cm4ge051bWJlcn0gICBBbmdsZSBpbiByYWRpYW5zXG4gICAqL1xuICBmdW5jdGlvbiB0b3JhZChkKSB7XG4gICAgcmV0dXJuIGQgKiBNYXRoLlBJIC8gMTgwLjA7XG4gIH1cblxuICAvKipcbiAgICogQ29udmVydCByYWRpYW5zIHRvIGRlZ3JlZXNcbiAgICogQHBhcmFtICB7TnVtYmVyfSByIEFuZ2xlIGluIHJhZGlhbnNcbiAgICogQHJldHVybiB7TnVtYmVyfSAgIEFuZ2xlIGluIGRlZ3JlZXNcbiAgICovXG4gIGZ1bmN0aW9uIHRvZGVnKHIpIHtcbiAgICByZXR1cm4gciAqIDE4MC4wIC8gTWF0aC5QSTtcbiAgfVxuXG4gIGZ1bmN0aW9uIGRzaW4oZCkge1xuICAgIHJldHVybiBNYXRoLnNpbih0b3JhZChkKSk7XG4gIH1cblxuICBmdW5jdGlvbiBkY29zKGQpIHtcbiAgICByZXR1cm4gTWF0aC5jb3ModG9yYWQoZCkpO1xuICB9XG5cbiAgLyoqXG4gICAqIFNvbHZlIHRoZSBlcXVhdGlvbiBvZiBLZXBsZXIuXG4gICAqL1xuICBmdW5jdGlvbiBrZXBsZXIobSwgZWNjKSB7XG4gICAgdmFyIGVwc2lsb24gPSAxZS02O1xuXG4gICAgbSA9IHRvcmFkKG0pO1xuICAgIHZhciBlID0gbTtcbiAgICB3aGlsZSgxKSB7XG4gICAgICB2YXIgZGVsdGEgPSBlIC0gZWNjICogTWF0aC5zaW4oZSkgLSBtO1xuICAgICAgZSA9IGUgLSBkZWx0YSAvICgxLjAgLSBlY2MgKiBNYXRoLmNvcyhlKSk7XG5cbiAgICAgIGlmIChNYXRoLmFicyhkZWx0YSkgPD0gZXBzaWxvbikge1xuICAgICAgICBicmVhaztcbiAgICAgIH1cbiAgICB9XG5cbiAgICByZXR1cm4gZTtcbiAgfVxuXG4gIC8qKlxuICAgKiBGaW5kcyB0aGUgcGhhc2UgaW5mb3JtYXRpb24gZm9yIHNwZWNpZmljIGRhdGUuXG4gICAqIEBwYXJhbSAge0RhdGV9IHBoYXNlX2RhdGUgRGF0ZSB0byBnZXQgcGhhc2UgaW5mb3JtYXRpb24gb2YuXG4gICAqIEByZXR1cm4ge09iamVjdH0gICAgICAgICAgUGhhc2UgZGF0YVxuICAgKi9cbiAgZnVuY3Rpb24gcGhhc2UocGhhc2VfZGF0ZSkge1xuICAgIGlmKCFwaGFzZV9kYXRlKSB7XG4gICAgICBwaGFzZV9kYXRlID0gKG5ldyBEYXRlKCkpLmdldEp1bGlhbigpO1xuICAgIH1cbiAgICBlbHNlIHtcbiAgICAgIHBoYXNlX2RhdGUgPSBwaGFzZV9kYXRlLmdldEp1bGlhbigpO1xuICAgIH1cblxuICAgIHZhciBkYXkgPSBwaGFzZV9kYXRlIC0gYy5lcG9jaDtcblxuICAgIC8vIE1lYW4gYW5vbWFseSBvZiB0aGUgU3VuXG4gICAgdmFyIE4gPSBmaXhhbmdsZSgoMzYwLzM2NS4yNDIyKSAqIGRheSk7XG4gICAgLy9Db252ZXJ0IGZyb20gcGVyaWdlZSBjb29yZGluYXRlcyB0byBlcG9jaCAxOTgwXG4gICAgdmFyIE0gPSBmaXhhbmdsZShOICsgYy5lY2xpcHRpY19sb25naXR1ZGVfZXBvY2ggLSBjLmVjbGlwdGljX2xvbmdpdHVkZV9wZXJpZ2VlKTtcblxuICAgIC8vIFNvbHZlIEtlcGxlcidzIGVxdWF0aW9uXG4gICAgdmFyIEVjID0ga2VwbGVyKE0sIGMuZWNjZW50cmljaXR5KTtcbiAgICBFYyA9IE1hdGguc3FydCgoMSArIGMuZWNjZW50cmljaXR5KSAvICgxIC0gYy5lY2NlbnRyaWNpdHkpKSAqIE1hdGgudGFuKEVjLzIuMCk7XG4gICAgLy8gVHJ1ZSBhbm9tYWx5XG4gICAgRWMgPSAyICogdG9kZWcoTWF0aC5hdGFuKEVjKSk7XG4gICAgLy8gU3VucydzIGdlb21ldHJpYyBlY2xpcHRpYyBsb25ndWl0dWRlXG4gICAgdmFyIGxhbWJkYV9zdW4gPSBmaXhhbmdsZShFYyArIGMuZWNsaXB0aWNfbG9uZ2l0dWRlX3BlcmlnZWUpO1xuXG4gICAgLy8gT3JiaXRhbCBkaXN0YW5jZSBmYWN0b3JcbiAgICB2YXIgRiA9ICgoMSArIGMuZWNjZW50cmljaXR5ICogTWF0aC5jb3ModG9yYWQoRWMpKSkgLyAoMSAtIE1hdGgucG93KGMuZWNjZW50cmljaXR5LCAyKSkpO1xuXG4gICAgLy8gRGlzdGFuY2UgdG8gU3VuIGluIGttXG4gICAgdmFyIHN1bl9kaXN0ID0gYy5zdW5fc21heGlzIC8gRjtcbiAgICB2YXIgc3VuX2FuZ3VsYXJfZGlhbWV0ZXIgPSBGICogYy5zdW5fYW5ndWxhcl9zaXplX3NtYXhpcztcblxuICAgIC8vIENhbGN1bGF0aW9uIG9mIHRoZSBNb29uJ3MgcG9zaXRpb25cblxuICAgIC8vIE1vb24ncyBtZWFuIGxvbmdpdHVkZVxuICAgIHZhciBtb29uX2xvbmdpdHVkZSA9IGZpeGFuZ2xlKDEzLjE3NjM5NjYgKiBkYXkgKyBjLm1vb25fbWVhbl9sb25naXR1ZGVfZXBvY2gpO1xuXG4gICAgLy8gTW9vbidzIG1lYW4gYW5vbWFseVxuICAgIHZhciBNTSA9IGZpeGFuZ2xlKG1vb25fbG9uZ2l0dWRlIC0gMC4xMTE0MDQxICogZGF5IC0gYy5tb29uX21lYW5fcGVyaWdlZV9lcG9jaCk7XG5cbiAgICAvLyBNb29uJ3MgYXNjZW5kaW5nIG5vZGUgbWVhbiBsb25naXR1ZGVcbiAgICAvLyBNTiA9IGZpeGFuZ2xlKGMubm9kZV9tZWFuX2xvbmdpdHVkZV9lcG9jaCAtIDAuMDUyOTUzOSAqIGRheSlcblxuICAgIHZhciBldmVjdGlvbiA9IDEuMjczOSAqIE1hdGguc2luKHRvcmFkKDIqKG1vb25fbG9uZ2l0dWRlIC0gbGFtYmRhX3N1bikgLSBNTSkpO1xuXG4gICAgLy8gQW5udWFsIGVxdWF0aW9uXG4gICAgdmFyIGFubnVhbF9lcSA9IDAuMTg1OCAqIE1hdGguc2luKHRvcmFkKE0pKTtcblxuICAgIC8vIENvcnJlY3Rpb24gdGVybVxuICAgIHZhciBBMyA9IDAuMzcgKiBNYXRoLnNpbih0b3JhZChNKSk7XG5cbiAgICB2YXIgTW1QID0gTU0gKyBldmVjdGlvbiAtIGFubnVhbF9lcSAtIEEzO1xuXG4gICAgLy8gQ29ycmVjdGlvbiBmb3IgdGhlIGVxdWF0aW9uIG9mIHRoZSBjZW50cmVcbiAgICB2YXIgbUVjID0gNi4yODg2ICogTWF0aC5zaW4odG9yYWQoTW1QKSk7XG5cbiAgICAvLyBBbm90aGVyIGNvcnJlY3Rpb24gdGVybVxuICAgIHZhciBBNCA9IDAuMjE0ICogTWF0aC5zaW4odG9yYWQoMiAqIE1tUCkpO1xuXG4gICAgLy8gQ29ycmVjdGVkIGxvbmdpdHVkZVxuICAgIHZhciBsUCA9IG1vb25fbG9uZ2l0dWRlICsgZXZlY3Rpb24gKyBtRWMgLSBhbm51YWxfZXEgKyBBNDtcblxuICAgIC8vIFZhcmlhdGlvblxuICAgIHZhciB2YXJpYXRpb24gPSAwLjY1ODMgKiBNYXRoLnNpbih0b3JhZCgyKihsUCAtIGxhbWJkYV9zdW4pKSk7XG5cbiAgICAvLyBUcnVlIGxvbmdpdHVkZVxuICAgIHZhciBsUFAgPSBsUCArIHZhcmlhdGlvbjtcblxuICAgIC8vIENhbGN1bGF0aW9uIG9mIHRoZSBwaGFzZSBvZiB0aGUgTW9vblxuXG4gICAgLy8gQWdlIG9mIHRoZSBNb29uLCBpbiBkZWdyZWVzXG4gICAgdmFyIG1vb25fYWdlID0gbFBQIC0gbGFtYmRhX3N1bjtcblxuICAgIC8vIFBoYXNlIG9mIHRoZSBNb29uXG4gICAgdmFyIG1vb25fcGhhc2UgPSAoMSAtIE1hdGguY29zKHRvcmFkKG1vb25fYWdlKSkpIC8gMi4wO1xuXG4gICAgLy8gQ2FsY3VsYXRlIGRpc3RhbmNlIG9mIE1vb24gZnJvbSB0aGUgY2VudHJlIG9mIHRoZSBFYXJ0aFxuICAgIHZhciBtb29uX2Rpc3QgPSAoYy5tb29uX3NtYXhpcyAqICgxIC0gTWF0aC5wb3coYy5tb29uX2VjY2VudHJpY2l0eSwyKSkpIC8gKDEgKyBjLm1vb25fZWNjZW50cmljaXR5ICogTWF0aC5jb3ModG9yYWQoTW1QICsgbUVjKSkpO1xuXG4gICAgLy8gQ2FsY3VsYXRlIE1vb24ncyBhbmd1bGFyIGRpYW1ldGVyXG4gICAgdmFyIG1vb25fZGlhbV9mcmFjID0gbW9vbl9kaXN0IC8gYy5tb29uX3NtYXhpcztcbiAgICB2YXIgbW9vbl9hbmd1bGFyX2RpYW1ldGVyID0gYy5tb29uX2FuZ3VsYXJfc2l6ZSAvIG1vb25fZGlhbV9mcmFjO1xuXG4gICAgLy8gQ2FsY3VsYXRlIE1vb24ncyBwYXJhbGxheCAodW51c2VkPylcbiAgICAvLyBtb29uX3BhcmFsbGF4ID0gYy5tb29uX3BhcmFsbGF4IC8gbW9vbl9kaWFtX2ZyYWNcblxuICAgIHZhciByZXMgPSB7XG4gICAgICAncGhhc2UnOiBmaXhhbmdsZShtb29uX2FnZSkgLyAzNjAuMCxcbiAgICAgICdpbGx1bWluYXRlZCc6IG1vb25fcGhhc2UsXG4gICAgICAnYWdlJzogYy5zeW5vZGljX21vbnRoICogZml4YW5nbGUobW9vbl9hZ2UpIC8gMzYwLjAsXG4gICAgICAnZGlzdGFuY2UnOiBtb29uX2Rpc3QsXG4gICAgICAnYW5ndWxhcl9kaWFtZXRlcic6IG1vb25fYW5ndWxhcl9kaWFtZXRlcixcbiAgICAgICdzdW5fZGlzdGFuY2UnOiBzdW5fZGlzdCxcbiAgICAgICdzdW5fYW5ndWxhcl9kaWFtZXRlcic6IHN1bl9hbmd1bGFyX2RpYW1ldGVyXG4gICAgfTtcblxuICAgIHJldHVybiByZXM7XG4gIH1cblxuICAvKipcbiAgICogRmluZCB0aW1lIG9mIHBoYXNlcyBvZiB0aGUgbW9vbiB3aGljaCBzdXJyb3VuZCB0aGUgY3VycmVudCBkYXRlLlxuICAgKiBGaXZlIHBoYXNlcyBhcmUgZm91bmQsIHN0YXJ0aW5nIGFuZCBlbmRpbmcgd2l0aCB0aGUgbmV3IG1vb25zXG4gICAqIHdoaWNoIGJvdW5kIHRoZSBjdXJyZW50IGx1bmF0aW9uLlxuICAgKiBAcGFyYW0gIHtEYXRlfSBzZGF0ZSBEYXRlIHRvIHN0YXJ0IGh1bnRpbmcgZnJvbSAoZGVmYXVsdHMgdG8gY3VycmVudCBkYXRlKVxuICAgKiBAcmV0dXJuIHtPYmplY3R9ICAgICBPYmplY3QgY29udGFpbmluZyByZWNlbnQgcGFzdCBhbmQgZnV0dXJlIHBoYXNlc1xuICAgKi9cbiAgZnVuY3Rpb24gcGhhc2VfaHVudChzZGF0ZSkge1xuICAgIGlmKCFzZGF0ZSkge1xuICAgICAgc2RhdGUgPSBuZXcgRGF0ZSgpO1xuICAgIH1cblxuICAgIHZhciBhZGF0ZSA9IG5ldyBEYXRlKHNkYXRlLnZhbHVlT2YoKSk7IC8vIHRvZGF5IVxuICAgIHZhciB4ID0gNDU7IC8vIGdvIGJhY2sgNDUgZGF5cyFcbiAgICBhZGF0ZS5zZXREYXRlKGFkYXRlLmdldERhdGUoKSAtIHgpO1xuXG4gICAgdmFyIGsxID0gTWF0aC5mbG9vcigoYWRhdGUuZ2V0RnVsbFllYXIoKSArICgoYWRhdGUuZ2V0TW9udGgoKSkgKiAoMS4wLzEyLjApKSAtIDE5MDApICogMTIuMzY4NSk7XG4gICAgdmFyIG50MSA9IG1lYW5waGFzZShhZGF0ZSwgazEpO1xuICAgIGFkYXRlID0gbnQxO1xuXG4gICAgc2RhdGUgPSBzZGF0ZS5nZXRKdWxpYW4oKTtcbiAgICB2YXIgazI7XG4gICAgd2hpbGUoMSkge1xuICAgICAgYWRhdGUgPSBhZGF0ZSArIGMuc3lub2RpY19tb250aDtcbiAgICAgIGsyID0gazEgKyAxO1xuICAgICAgdmFyIG50MiA9IG1lYW5waGFzZShhZGF0ZSwgazIpO1xuICAgICAgaWYobnQxIDw9IHNkYXRlICYmIHNkYXRlIDwgbnQyKSB7XG4gICAgICAgIGJyZWFrO1xuICAgICAgfVxuICAgICAgbnQxID0gbnQyO1xuICAgICAgazEgPSBrMjtcbiAgICB9XG4gICAgdmFyIGtzID0gW2sxLCBrMSwgazEsIGsxLCBrMl07XG4gICAgdmFyIHRwaGFzZXMgPSBbTkVXLCBGSVJTVCwgRlVMTCwgTEFTVCwgTkVXXTtcbiAgICB2YXIgcGhhc2VfbmFtZXMgPSBbJ25ld19kYXRlJywgJ3ExX2RhdGUnLCAnZnVsbF9kYXRlJywgJ3EzX2RhdGUnLCAnbmV4dG5ld19kYXRlJ107XG4gICAgdmFyIHBoYXNlcyA9IHt9O1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBrcy5sZW5ndGg7IGkrKykge1xuICAgICAgcGhhc2VzW3BoYXNlX25hbWVzW2ldXSA9IHRydWVwaGFzZShrc1tpXSwgdHBoYXNlc1tpXSk7XG4gICAgfVxuXG4gICAgcmV0dXJuIHBoYXNlcztcbiAgfVxuXG4gIC8qKlxuICAgKiBHaXZlbiBhIEsgdmFsdWUgdXNlZCB0byBkZXRlcm1pbmUgdGhlIG1lYW4gcGhhc2Ugb2YgdGhlIG5ld1xuICAgKiBtb29uLCBhbmQgYSBwaGFzZSBzZWxlY3RvciAoMC4wLCAwLjI1LCAwLjUsIDAuNzUpLCBvYnRhaW4gdGhlXG4gICAqIHRydWUsIGNvcnJlY3RlZCBwaGFzZSB0aW1lLlxuICAgKiBAcGFyYW0gIHtbdHlwZV19IGsgICAgICBbZGVzY3JpcHRpb25dXG4gICAqIEBwYXJhbSAge1t0eXBlXX0gdHBoYXNlIFtkZXNjcmlwdGlvbl1cbiAgICogQHJldHVybiB7W3R5cGVdfSAgICAgICAgW2Rlc2NyaXB0aW9uXVxuICAgKi9cbiAgZnVuY3Rpb24gdHJ1ZXBoYXNlKGssIHRwaGFzZSkge1xuXG4gICAgdmFyIGFwY29yID0gZmFsc2U7XG5cbiAgICAvLyBhZGQgcGhhc2UgdG8gbmV3IG1vb24gdGltZVxuICAgIGsgPSBrICsgdHBoYXNlO1xuICAgIC8vIFRpbWUgaW4gSnVsaWFuIGNlbnR1cmllcyBmcm9tIDE5MDAgSmFudWFyeSAwLjVcbiAgICB2YXIgdCA9IGsgLyAxMjM2Ljg1O1xuXG4gICAgdmFyIHQyID0gdCAqIHQ7XG4gICAgdmFyIHQzID0gdDIgKiB0O1xuXG4gICAgLy8gTWVhbiB0aW1lIG9mIHBoYXNlXG4gICAgdmFyIHB0ID0gKFxuICAgICAgMjQxNTAyMC43NTkzMyArIGMuc3lub2RpY19tb250aCAqIGsgKyAwLjAwMDExNzggKiB0MiAtXG4gICAgICAwLjAwMDAwMDE1NSAqIHQzICsgMC4wMDAzMyAqIGRzaW4oMTY2LjU2ICsgMTMyLjg3ICogdCAtXG4gICAgICAwLjAwOTE3MyAqIHQyKVxuICAgICk7XG5cbiAgICAvLyBTdW4ncyBtZWFuIGFub21hbHlcbiAgICB2YXIgbSA9IDM1OS4yMjQyICsgMjkuMTA1MzU2MDggKiBrIC0gMC4wMDAwMzMzICogdDIgLSAwLjAwMDAwMzQ3ICogdDM7XG5cbiAgICAvLyBNb29uJ3MgbWVhbiBhbm9tYWx5XG4gICAgdmFyIG1wcmltZSA9IDMwNi4wMjUzICsgMzg1LjgxNjkxODA2ICogayArIDAuMDEwNzMwNiAqIHQyICsgMC4wMDAwMTIzNiAqIHQzO1xuXG4gICAgLy8gTW9vbidzIGFyZ3VtZW50IG9mIGxhdGl0dWRlXG4gICAgdmFyIGYgPSAyMS4yOTY0ICsgMzkwLjY3MDUwNjQ2ICogayAtIDAuMDAxNjUyOCAqIHQyIC0gMC4wMDAwMDIzOSAqIHQzO1xuXG4gICAgaWYgKCh0cGhhc2UgPCAwLjAxKSB8fCAoTWF0aC5hYnModHBoYXNlIC0gMC41KSA8IDAuMDEpKSB7XG5cbiAgICAgIC8vIENvcnJlY3Rpb25zIGZvciBOZXcgYW5kIEZ1bGwgTW9vblxuICAgICAgcHQgPSBwdCArIChcbiAgICAgICAgKDAuMTczNCAtIDAuMDAwMzkzICogdCkgKiBkc2luKG0pICtcbiAgICAgICAgMC4wMDIxICogZHNpbigyICogbSkgLVxuICAgICAgICAwLjQwNjggKiBkc2luKG1wcmltZSkgK1xuICAgICAgICAwLjAxNjEgKiBkc2luKDIgKiBtcHJpbWUpIC1cbiAgICAgICAgMC4wMDA0ICogZHNpbigzICogbXByaW1lKSArXG4gICAgICAgIDAuMDEwNCAqIGRzaW4oMiAqIGYpIC1cbiAgICAgICAgMC4wMDUxICogZHNpbihtICsgbXByaW1lKSAtXG4gICAgICAgIDAuMDA3NCAqIGRzaW4obSAtIG1wcmltZSkgK1xuICAgICAgICAwLjAwMDQgKiBkc2luKDIgKiBmICsgbSkgLVxuICAgICAgICAwLjAwMDQgKiBkc2luKDIgKiBmIC0gbSkgLVxuICAgICAgICAwLjAwMDYgKiBkc2luKDIgKiBmICsgbXByaW1lKSArXG4gICAgICAgIDAuMDAxMCAqIGRzaW4oMiAqIGYgLSBtcHJpbWUpICtcbiAgICAgICAgMC4wMDA1ICogZHNpbihtICsgMiAqIG1wcmltZSlcbiAgICAgICk7XG5cbiAgICAgIGFwY29yID0gdHJ1ZTtcbiAgICB9XG4gICAgZWxzZSBpZiAoKE1hdGguYWJzKHRwaGFzZSAtIDAuMjUpIDwgMC4wMSkgfHwgKE1hdGguYWJzKHRwaGFzZSAtIDAuNzUpIDwgMC4wMSkpIHtcbiAgICAgICAgcHQgPSBwdCArIChcbiAgICAgICAgICAoMC4xNzIxIC0gMC4wMDA0ICogdCkgKiBkc2luKG0pICtcbiAgICAgICAgICAwLjAwMjEgKiBkc2luKDIgKiBtKSAtXG4gICAgICAgICAgMC42MjgwICogZHNpbihtcHJpbWUpICtcbiAgICAgICAgICAwLjAwODkgKiBkc2luKDIgKiBtcHJpbWUpIC1cbiAgICAgICAgICAwLjAwMDQgKiBkc2luKDMgKiBtcHJpbWUpICtcbiAgICAgICAgICAwLjAwNzkgKiBkc2luKDIgKiBmKSAtXG4gICAgICAgICAgMC4wMTE5ICogZHNpbihtICsgbXByaW1lKSAtXG4gICAgICAgICAgMC4wMDQ3ICogZHNpbihtIC0gbXByaW1lKSArXG4gICAgICAgICAgMC4wMDAzICogZHNpbigyICogZiArIG0pIC1cbiAgICAgICAgICAwLjAwMDQgKiBkc2luKDIgKiBmIC0gbSkgLVxuICAgICAgICAgIDAuMDAwNiAqIGRzaW4oMiAqIGYgKyBtcHJpbWUpICtcbiAgICAgICAgICAwLjAwMjEgKiBkc2luKDIgKiBmIC0gbXByaW1lKSArXG4gICAgICAgICAgMC4wMDAzICogZHNpbihtICsgMiAqIG1wcmltZSkgK1xuICAgICAgICAgIDAuMDAwNCAqIGRzaW4obSAtIDIgKiBtcHJpbWUpIC1cbiAgICAgICAgICAwLjAwMDMgKiBkc2luKDIgKiBtICsgbXByaW1lKVxuICAgICAgICApO1xuICAgICAgaWYgKHRwaGFzZSA8IDAuNSkge1xuICAgICAgICAgIC8vICBGaXJzdCBxdWFydGVyIGNvcnJlY3Rpb25cbiAgICAgICAgICBwdCA9IHB0ICsgMC4wMDI4IC0gMC4wMDA0ICogZGNvcyhtKSArIDAuMDAwMyAqIGRjb3MobXByaW1lKTtcbiAgICAgIH1cbiAgICAgIGVsc2Uge1xuICAgICAgICAgIC8vICBMYXN0IHF1YXJ0ZXIgY29ycmVjdGlvblxuICAgICAgICAgIHB0ID0gcHQgKyAtMC4wMDI4ICsgMC4wMDA0ICogZGNvcyhtKSAtIDAuMDAwMyAqIGRjb3MobXByaW1lKTtcbiAgICAgIH1cbiAgICAgIGFwY29yID0gdHJ1ZTtcbiAgICB9XG5cbiAgICBpZiAoIWFwY29yKSB7XG4gICAgICBjb25zb2xlLmxvZyhcIlRSVUVQSEFTRSBjYWxsZWQgd2l0aCBpbnZhbGlkIHBoYXNlIHNlbGVjdG9yIFwiLCB0cGhhc2UpO1xuICAgIH1cblxuICAgIHJldHVybiBwdC5KdWxpYW4yRGF0ZSh0cnVlKTtcbiAgfVxuXG4gIC8qKlxuICAgKiBDYWxjdWxhdGVzIHRpbWUgb2YgdGhlIG1lYW4gbmV3IE1vb24gZm9yIGEgZ2l2ZW4gYmFzZSBkYXRlLlxuICAgKiBUaGlzIGFyZ3VtZW50IEsgdG8gdGhpcyBmdW5jdGlvbiBpcyB0aGUgcHJlY29tcHV0ZWQgc3lub2RpYyBtb250aFxuICAgKiBpbmRleCwgZ2l2ZW4gYnk6XG4gICAqICAgSyA9ICh5ZWFyIC0gMTkwMCkgKiAxMi4zNjg1XG4gICAqIHdoZXJlIHllYXIgaXMgZXhwcmVzc2VkIGFzIGEgeWVhciBhbmQgZnJhY3Rpb25hbCB5ZWFyLlxuICAgKiBAcGFyYW0gIHtEYXRlfSBzZGF0ZSAgIFN0YXJ0IGRhdGVcbiAgICogQHBhcmFtICB7W3R5cGVdfSBrICAgICBbZGVzY3JpcHRpb25dXG4gICAqIEByZXR1cm4ge1t0eXBlXX0gICAgICAgW2Rlc2NyaXB0aW9uXVxuICAgKi9cbiAgZnVuY3Rpb24gbWVhbnBoYXNlKHNkYXRlLCBrKSB7XG5cbiAgICAvLyBUaW1lIGluIEp1bGlhbiBjZW50dXJpZXMgZnJvbSAxOTAwIEphbnVhcnkgMTIgbm9vblxuICAgIHZhciBkZWx0YV90ID0gKHNkYXRlIC0gKG5ldyBEYXRlKDE5MDAsMCwxLDEyKSkpIC8gKDEwMDAqNjAqNjAqMjQpO1xuICAgIHZhciB0ID0gZGVsdGFfdCAvIDM2NTI1O1xuXG4gICAgLy8gc3F1YXJlIGZvciBmcmVxdWVudCB1c2VcbiAgICB2YXIgdDIgPSB0ICogdDtcbiAgICAvLyBhbmQgY3ViZVxuICAgIHZhciB0MyA9IHQyICogdDtcblxuICAgIG50MSA9IChcbiAgICAgIDI0MTUwMjAuNzU5MzMgKyBjLnN5bm9kaWNfbW9udGggKiBrICsgMC4wMDAxMTc4ICogdDIgLVxuICAgICAgMC4wMDAwMDAxNTUgKiB0MyArIDAuMDAwMzMgKiBkc2luKDE2Ni41NiArIDEzMi44NyAqIHQgLVxuICAgICAgMC4wMDkxNzMgKiB0MilcbiAgICApO1xuXG4gICAgcmV0dXJuIG50MTtcbiAgfVxuXG4gIG1vZHVsZS5leHBvcnRzID0ge1xuICAgJ3BoYXNlX2h1bnQnOiBwaGFzZV9odW50LFxuICAgJ3BoYXNlJzogcGhhc2VcbiAgfTtcbn0pKCk7XG4iLCIvLyAoYykgRGVhbiBNY05hbWVlIDxkZWFuQGdtYWlsLmNvbT4sIDIwMTMuXG4vL1xuLy8gaHR0cHM6Ly9naXRodWIuY29tL2RlYW5tL29tZ2dpZlxuLy9cbi8vIFBlcm1pc3Npb24gaXMgaGVyZWJ5IGdyYW50ZWQsIGZyZWUgb2YgY2hhcmdlLCB0byBhbnkgcGVyc29uIG9idGFpbmluZyBhIGNvcHlcbi8vIG9mIHRoaXMgc29mdHdhcmUgYW5kIGFzc29jaWF0ZWQgZG9jdW1lbnRhdGlvbiBmaWxlcyAodGhlIFwiU29mdHdhcmVcIiksIHRvXG4vLyBkZWFsIGluIHRoZSBTb2Z0d2FyZSB3aXRob3V0IHJlc3RyaWN0aW9uLCBpbmNsdWRpbmcgd2l0aG91dCBsaW1pdGF0aW9uIHRoZVxuLy8gcmlnaHRzIHRvIHVzZSwgY29weSwgbW9kaWZ5LCBtZXJnZSwgcHVibGlzaCwgZGlzdHJpYnV0ZSwgc3VibGljZW5zZSwgYW5kL29yXG4vLyBzZWxsIGNvcGllcyBvZiB0aGUgU29mdHdhcmUsIGFuZCB0byBwZXJtaXQgcGVyc29ucyB0byB3aG9tIHRoZSBTb2Z0d2FyZSBpc1xuLy8gZnVybmlzaGVkIHRvIGRvIHNvLCBzdWJqZWN0IHRvIHRoZSBmb2xsb3dpbmcgY29uZGl0aW9uczpcbi8vXG4vLyBUaGUgYWJvdmUgY29weXJpZ2h0IG5vdGljZSBhbmQgdGhpcyBwZXJtaXNzaW9uIG5vdGljZSBzaGFsbCBiZSBpbmNsdWRlZCBpblxuLy8gYWxsIGNvcGllcyBvciBzdWJzdGFudGlhbCBwb3J0aW9ucyBvZiB0aGUgU29mdHdhcmUuXG4vL1xuLy8gVEhFIFNPRlRXQVJFIElTIFBST1ZJREVEIFwiQVMgSVNcIiwgV0lUSE9VVCBXQVJSQU5UWSBPRiBBTlkgS0lORCwgRVhQUkVTUyBPUlxuLy8gSU1QTElFRCwgSU5DTFVESU5HIEJVVCBOT1QgTElNSVRFRCBUTyBUSEUgV0FSUkFOVElFUyBPRiBNRVJDSEFOVEFCSUxJVFksXG4vLyBGSVRORVNTIEZPUiBBIFBBUlRJQ1VMQVIgUFVSUE9TRSBBTkQgTk9OSU5GUklOR0VNRU5ULiBJTiBOTyBFVkVOVCBTSEFMTCBUSEVcbi8vIEFVVEhPUlMgT1IgQ09QWVJJR0hUIEhPTERFUlMgQkUgTElBQkxFIEZPUiBBTlkgQ0xBSU0sIERBTUFHRVMgT1IgT1RIRVJcbi8vIExJQUJJTElUWSwgV0hFVEhFUiBJTiBBTiBBQ1RJT04gT0YgQ09OVFJBQ1QsIFRPUlQgT1IgT1RIRVJXSVNFLCBBUklTSU5HXG4vLyBGUk9NLCBPVVQgT0YgT1IgSU4gQ09OTkVDVElPTiBXSVRIIFRIRSBTT0ZUV0FSRSBPUiBUSEUgVVNFIE9SIE9USEVSIERFQUxJTkdTXG4vLyBJTiBUSEUgU09GVFdBUkUuXG4vL1xuLy8gb21nZ2lmIGlzIGEgSmF2YVNjcmlwdCBpbXBsZW1lbnRhdGlvbiBvZiBhIEdJRiA4OWEgZW5jb2RlciBhbmQgZGVjb2Rlcixcbi8vIGluY2x1ZGluZyBhbmltYXRpb24gYW5kIGNvbXByZXNzaW9uLiAgSXQgZG9lcyBub3QgcmVseSBvbiBhbnkgc3BlY2lmaWNcbi8vIHVuZGVybHlpbmcgc3lzdGVtLCBzbyBzaG91bGQgcnVuIGluIHRoZSBicm93c2VyLCBOb2RlLCBvciBQbGFzay5cblxuZnVuY3Rpb24gR2lmV3JpdGVyKGJ1Ziwgd2lkdGgsIGhlaWdodCwgZ29wdHMpIHtcbiAgdmFyIHAgPSAwO1xuXG4gIHZhciBnb3B0cyA9IGdvcHRzID09PSB1bmRlZmluZWQgPyB7IH0gOiBnb3B0cztcbiAgdmFyIGxvb3BfY291bnQgPSBnb3B0cy5sb29wID09PSB1bmRlZmluZWQgPyBudWxsIDogZ29wdHMubG9vcDtcbiAgdmFyIGdsb2JhbF9wYWxldHRlID0gZ29wdHMucGFsZXR0ZSA9PT0gdW5kZWZpbmVkID8gbnVsbCA6IGdvcHRzLnBhbGV0dGU7XG5cbiAgaWYgKHdpZHRoIDw9IDAgfHwgaGVpZ2h0IDw9IDAgfHwgd2lkdGggPiA2NTUzNSB8fCBoZWlnaHQgPiA2NTUzNSlcbiAgICB0aHJvdyBcIldpZHRoL0hlaWdodCBpbnZhbGlkLlwiXG5cbiAgZnVuY3Rpb24gY2hlY2tfcGFsZXR0ZV9hbmRfbnVtX2NvbG9ycyhwYWxldHRlKSB7XG4gICAgdmFyIG51bV9jb2xvcnMgPSBwYWxldHRlLmxlbmd0aDtcbiAgICBpZiAobnVtX2NvbG9ycyA8IDIgfHwgbnVtX2NvbG9ycyA+IDI1NiB8fCAgbnVtX2NvbG9ycyAmIChudW1fY29sb3JzLTEpKVxuICAgICAgdGhyb3cgXCJJbnZhbGlkIGNvZGUvY29sb3IgbGVuZ3RoLCBtdXN0IGJlIHBvd2VyIG9mIDIgYW5kIDIgLi4gMjU2LlwiO1xuICAgIHJldHVybiBudW1fY29sb3JzO1xuICB9XG5cbiAgLy8gLSBIZWFkZXIuXG4gIGJ1ZltwKytdID0gMHg0NzsgYnVmW3ArK10gPSAweDQ5OyBidWZbcCsrXSA9IDB4NDY7ICAvLyBHSUZcbiAgYnVmW3ArK10gPSAweDM4OyBidWZbcCsrXSA9IDB4Mzk7IGJ1ZltwKytdID0gMHg2MTsgIC8vIDg5YVxuXG4gIC8vIEhhbmRsaW5nIG9mIEdsb2JhbCBDb2xvciBUYWJsZSAocGFsZXR0ZSkgYW5kIGJhY2tncm91bmQgaW5kZXguXG4gIHZhciBncF9udW1fY29sb3JzX3BvdzIgPSAwO1xuICB2YXIgYmFja2dyb3VuZCA9IDA7XG4gIGlmIChnbG9iYWxfcGFsZXR0ZSAhPT0gbnVsbCkge1xuICAgIHZhciBncF9udW1fY29sb3JzID0gY2hlY2tfcGFsZXR0ZV9hbmRfbnVtX2NvbG9ycyhnbG9iYWxfcGFsZXR0ZSk7XG4gICAgd2hpbGUgKGdwX251bV9jb2xvcnMgPj49IDEpICsrZ3BfbnVtX2NvbG9yc19wb3cyO1xuICAgIGdwX251bV9jb2xvcnMgPSAxIDw8IGdwX251bV9jb2xvcnNfcG93MjtcbiAgICAtLWdwX251bV9jb2xvcnNfcG93MjtcbiAgICBpZiAoZ29wdHMuYmFja2dyb3VuZCAhPT0gdW5kZWZpbmVkKSB7XG4gICAgICBiYWNrZ3JvdW5kID0gZ29wdHMuYmFja2dyb3VuZDtcbiAgICAgIGlmIChiYWNrZ3JvdW5kID49IGdwX251bV9jb2xvcnMpIHRocm93IFwiQmFja2dyb3VuZCBpbmRleCBvdXQgb2YgcmFuZ2UuXCI7XG4gICAgICAvLyBUaGUgR0lGIHNwZWMgc3RhdGVzIHRoYXQgYSBiYWNrZ3JvdW5kIGluZGV4IG9mIDAgc2hvdWxkIGJlIGlnbm9yZWQsIHNvXG4gICAgICAvLyB0aGlzIGlzIHByb2JhYmx5IGEgbWlzdGFrZSBhbmQgeW91IHJlYWxseSB3YW50IHRvIHNldCBpdCB0byBhbm90aGVyXG4gICAgICAvLyBzbG90IGluIHRoZSBwYWxldHRlLiAgQnV0IGFjdHVhbGx5IGluIHRoZSBlbmQgbW9zdCBicm93c2VycywgZXRjIGVuZFxuICAgICAgLy8gdXAgaWdub3JpbmcgdGhpcyBhbG1vc3QgY29tcGxldGVseSAoaW5jbHVkaW5nIGZvciBkaXNwb3NlIGJhY2tncm91bmQpLlxuICAgICAgaWYgKGJhY2tncm91bmQgPT09IDApXG4gICAgICAgIHRocm93IFwiQmFja2dyb3VuZCBpbmRleCBleHBsaWNpdGx5IHBhc3NlZCBhcyAwLlwiO1xuICAgIH1cbiAgfVxuXG4gIC8vIC0gTG9naWNhbCBTY3JlZW4gRGVzY3JpcHRvci5cbiAgLy8gTk9URShkZWFubSk6IHcvaCBhcHBhcmVudGx5IGlnbm9yZWQgYnkgaW1wbGVtZW50YXRpb25zLCBidXQgc2V0IGFueXdheS5cbiAgYnVmW3ArK10gPSB3aWR0aCAmIDB4ZmY7IGJ1ZltwKytdID0gd2lkdGggPj4gOCAmIDB4ZmY7XG4gIGJ1ZltwKytdID0gaGVpZ2h0ICYgMHhmZjsgYnVmW3ArK10gPSBoZWlnaHQgPj4gOCAmIDB4ZmY7XG4gIC8vIE5PVEU6IEluZGljYXRlcyAwLWJwcCBvcmlnaW5hbCBjb2xvciByZXNvbHV0aW9uICh1bnVzZWQ/KS5cbiAgYnVmW3ArK10gPSAoZ2xvYmFsX3BhbGV0dGUgIT09IG51bGwgPyAweDgwIDogMCkgfCAgLy8gR2xvYmFsIENvbG9yIFRhYmxlIEZsYWcuXG4gICAgICAgICAgICAgZ3BfbnVtX2NvbG9yc19wb3cyOyAgLy8gTk9URTogTm8gc29ydCBmbGFnICh1bnVzZWQ/KS5cbiAgYnVmW3ArK10gPSBiYWNrZ3JvdW5kOyAgLy8gQmFja2dyb3VuZCBDb2xvciBJbmRleC5cbiAgYnVmW3ArK10gPSAwOyAgLy8gUGl4ZWwgYXNwZWN0IHJhdGlvICh1bnVzZWQ/KS5cblxuICAvLyAtIEdsb2JhbCBDb2xvciBUYWJsZVxuICBpZiAoZ2xvYmFsX3BhbGV0dGUgIT09IG51bGwpIHtcbiAgICBmb3IgKHZhciBpID0gMCwgaWwgPSBnbG9iYWxfcGFsZXR0ZS5sZW5ndGg7IGkgPCBpbDsgKytpKSB7XG4gICAgICB2YXIgcmdiID0gZ2xvYmFsX3BhbGV0dGVbaV07XG4gICAgICBidWZbcCsrXSA9IHJnYiA+PiAxNiAmIDB4ZmY7XG4gICAgICBidWZbcCsrXSA9IHJnYiA+PiA4ICYgMHhmZjtcbiAgICAgIGJ1ZltwKytdID0gcmdiICYgMHhmZjtcbiAgICB9XG4gIH1cblxuICBpZiAobG9vcF9jb3VudCAhPT0gbnVsbCkgeyAgLy8gTmV0c2NhcGUgYmxvY2sgZm9yIGxvb3BpbmcuXG4gICAgaWYgKGxvb3BfY291bnQgPCAwIHx8IGxvb3BfY291bnQgPiA2NTUzNSlcbiAgICAgIHRocm93IFwiTG9vcCBjb3VudCBpbnZhbGlkLlwiXG4gICAgLy8gRXh0ZW5zaW9uIGNvZGUsIGxhYmVsLCBhbmQgbGVuZ3RoLlxuICAgIGJ1ZltwKytdID0gMHgyMTsgYnVmW3ArK10gPSAweGZmOyBidWZbcCsrXSA9IDB4MGI7XG4gICAgLy8gTkVUU0NBUEUyLjBcbiAgICBidWZbcCsrXSA9IDB4NGU7IGJ1ZltwKytdID0gMHg0NTsgYnVmW3ArK10gPSAweDU0OyBidWZbcCsrXSA9IDB4NTM7XG4gICAgYnVmW3ArK10gPSAweDQzOyBidWZbcCsrXSA9IDB4NDE7IGJ1ZltwKytdID0gMHg1MDsgYnVmW3ArK10gPSAweDQ1O1xuICAgIGJ1ZltwKytdID0gMHgzMjsgYnVmW3ArK10gPSAweDJlOyBidWZbcCsrXSA9IDB4MzA7XG4gICAgLy8gU3ViLWJsb2NrXG4gICAgYnVmW3ArK10gPSAweDAzOyBidWZbcCsrXSA9IDB4MDE7XG4gICAgYnVmW3ArK10gPSBsb29wX2NvdW50ICYgMHhmZjsgYnVmW3ArK10gPSBsb29wX2NvdW50ID4+IDggJiAweGZmO1xuICAgIGJ1ZltwKytdID0gMHgwMDsgIC8vIFRlcm1pbmF0b3IuXG4gIH1cblxuXG4gIHZhciBlbmRlZCA9IGZhbHNlO1xuXG4gIHRoaXMuYWRkRnJhbWUgPSBmdW5jdGlvbih4LCB5LCB3LCBoLCBpbmRleGVkX3BpeGVscywgb3B0cykge1xuICAgIGlmIChlbmRlZCA9PT0gdHJ1ZSkgeyAtLXA7IGVuZGVkID0gZmFsc2U7IH0gIC8vIFVuLWVuZC5cblxuICAgIG9wdHMgPSBvcHRzID09PSB1bmRlZmluZWQgPyB7IH0gOiBvcHRzO1xuXG4gICAgLy8gVE9ETyhkZWFubSk6IEJvdW5kcyBjaGVjayB4LCB5LiAgRG8gdGhleSBuZWVkIHRvIGJlIHdpdGhpbiB0aGUgdmlydHVhbFxuICAgIC8vIGNhbnZhcyB3aWR0aC9oZWlnaHQsIEkgaW1hZ2luZT9cbiAgICBpZiAoeCA8IDAgfHwgeSA8IDAgfHwgeCA+IDY1NTM1IHx8IHkgPiA2NTUzNSlcbiAgICAgIHRocm93IFwieC95IGludmFsaWQuXCJcblxuICAgIGlmICh3IDw9IDAgfHwgaCA8PSAwIHx8IHcgPiA2NTUzNSB8fCBoID4gNjU1MzUpXG4gICAgICB0aHJvdyBcIldpZHRoL0hlaWdodCBpbnZhbGlkLlwiXG5cbiAgICBpZiAoaW5kZXhlZF9waXhlbHMubGVuZ3RoIDwgdyAqIGgpXG4gICAgICB0aHJvdyBcIk5vdCBlbm91Z2ggcGl4ZWxzIGZvciB0aGUgZnJhbWUgc2l6ZS5cIjtcblxuICAgIHZhciB1c2luZ19sb2NhbF9wYWxldHRlID0gdHJ1ZTtcbiAgICB2YXIgcGFsZXR0ZSA9IG9wdHMucGFsZXR0ZTtcbiAgICBpZiAocGFsZXR0ZSA9PT0gdW5kZWZpbmVkIHx8IHBhbGV0dGUgPT09IG51bGwpIHtcbiAgICAgIHVzaW5nX2xvY2FsX3BhbGV0dGUgPSBmYWxzZTtcbiAgICAgIHBhbGV0dGUgPSBnbG9iYWxfcGFsZXR0ZTtcbiAgICB9XG5cbiAgICBpZiAocGFsZXR0ZSA9PT0gdW5kZWZpbmVkIHx8IHBhbGV0dGUgPT09IG51bGwpXG4gICAgICB0aHJvdyBcIk11c3Qgc3VwcGx5IGVpdGhlciBhIGxvY2FsIG9yIGdsb2JhbCBwYWxldHRlLlwiO1xuXG4gICAgdmFyIG51bV9jb2xvcnMgPSBjaGVja19wYWxldHRlX2FuZF9udW1fY29sb3JzKHBhbGV0dGUpO1xuXG4gICAgLy8gQ29tcHV0ZSB0aGUgbWluX2NvZGVfc2l6ZSAocG93ZXIgb2YgMiksIGRlc3Ryb3lpbmcgbnVtX2NvbG9ycy5cbiAgICB2YXIgbWluX2NvZGVfc2l6ZSA9IDA7XG4gICAgd2hpbGUgKG51bV9jb2xvcnMgPj49IDEpICsrbWluX2NvZGVfc2l6ZTtcbiAgICBudW1fY29sb3JzID0gMSA8PCBtaW5fY29kZV9zaXplOyAgLy8gTm93IHdlIGNhbiBlYXNpbHkgZ2V0IGl0IGJhY2suXG5cbiAgICB2YXIgZGVsYXkgPSBvcHRzLmRlbGF5ID09PSB1bmRlZmluZWQgPyAwIDogb3B0cy5kZWxheTtcblxuICAgIC8vIEZyb20gdGhlIHNwZWM6XG4gICAgLy8gICAgIDAgLSAgIE5vIGRpc3Bvc2FsIHNwZWNpZmllZC4gVGhlIGRlY29kZXIgaXNcbiAgICAvLyAgICAgICAgICAgbm90IHJlcXVpcmVkIHRvIHRha2UgYW55IGFjdGlvbi5cbiAgICAvLyAgICAgMSAtICAgRG8gbm90IGRpc3Bvc2UuIFRoZSBncmFwaGljIGlzIHRvIGJlIGxlZnRcbiAgICAvLyAgICAgICAgICAgaW4gcGxhY2UuXG4gICAgLy8gICAgIDIgLSAgIFJlc3RvcmUgdG8gYmFja2dyb3VuZCBjb2xvci4gVGhlIGFyZWEgdXNlZCBieSB0aGVcbiAgICAvLyAgICAgICAgICAgZ3JhcGhpYyBtdXN0IGJlIHJlc3RvcmVkIHRvIHRoZSBiYWNrZ3JvdW5kIGNvbG9yLlxuICAgIC8vICAgICAzIC0gICBSZXN0b3JlIHRvIHByZXZpb3VzLiBUaGUgZGVjb2RlciBpcyByZXF1aXJlZCB0b1xuICAgIC8vICAgICAgICAgICByZXN0b3JlIHRoZSBhcmVhIG92ZXJ3cml0dGVuIGJ5IHRoZSBncmFwaGljIHdpdGhcbiAgICAvLyAgICAgICAgICAgd2hhdCB3YXMgdGhlcmUgcHJpb3IgdG8gcmVuZGVyaW5nIHRoZSBncmFwaGljLlxuICAgIC8vICA0LTcgLSAgICBUbyBiZSBkZWZpbmVkLlxuICAgIC8vIE5PVEUoZGVhbm0pOiBEaXNwb3NlIGJhY2tncm91bmQgZG9lc24ndCByZWFsbHkgd29yaywgYXBwYXJlbnRseSBtb3N0XG4gICAgLy8gYnJvd3NlcnMgaWdub3JlIHRoZSBiYWNrZ3JvdW5kIHBhbGV0dGUgaW5kZXggYW5kIGNsZWFyIHRvIHRyYW5zcGFyZW5jeS5cbiAgICB2YXIgZGlzcG9zYWwgPSBvcHRzLmRpc3Bvc2FsID09PSB1bmRlZmluZWQgPyAwIDogb3B0cy5kaXNwb3NhbDtcbiAgICBpZiAoZGlzcG9zYWwgPCAwIHx8IGRpc3Bvc2FsID4gMykgIC8vIDQtNyBpcyByZXNlcnZlZC5cbiAgICAgIHRocm93IFwiRGlzcG9zYWwgb3V0IG9mIHJhbmdlLlwiO1xuXG4gICAgdmFyIHVzZV90cmFuc3BhcmVuY3kgPSBmYWxzZTtcbiAgICB2YXIgdHJhbnNwYXJlbnRfaW5kZXggPSAwO1xuICAgIGlmIChvcHRzLnRyYW5zcGFyZW50ICE9PSB1bmRlZmluZWQgJiYgb3B0cy50cmFuc3BhcmVudCAhPT0gbnVsbCkge1xuICAgICAgdXNlX3RyYW5zcGFyZW5jeSA9IHRydWU7XG4gICAgICB0cmFuc3BhcmVudF9pbmRleCA9IG9wdHMudHJhbnNwYXJlbnQ7XG4gICAgICBpZiAodHJhbnNwYXJlbnRfaW5kZXggPCAwIHx8IHRyYW5zcGFyZW50X2luZGV4ID49IG51bV9jb2xvcnMpXG4gICAgICAgIHRocm93IFwiVHJhbnNwYXJlbnQgY29sb3IgaW5kZXguXCI7XG4gICAgfVxuXG4gICAgaWYgKGRpc3Bvc2FsICE9PSAwIHx8IHVzZV90cmFuc3BhcmVuY3kgfHwgZGVsYXkgIT09IDApIHtcbiAgICAgIC8vIC0gR3JhcGhpY3MgQ29udHJvbCBFeHRlbnNpb25cbiAgICAgIGJ1ZltwKytdID0gMHgyMTsgYnVmW3ArK10gPSAweGY5OyAgLy8gRXh0ZW5zaW9uIC8gTGFiZWwuXG4gICAgICBidWZbcCsrXSA9IDQ7ICAvLyBCeXRlIHNpemUuXG5cbiAgICAgIGJ1ZltwKytdID0gZGlzcG9zYWwgPDwgMiB8ICh1c2VfdHJhbnNwYXJlbmN5ID09PSB0cnVlID8gMSA6IDApO1xuICAgICAgYnVmW3ArK10gPSBkZWxheSAmIDB4ZmY7IGJ1ZltwKytdID0gZGVsYXkgPj4gOCAmIDB4ZmY7XG4gICAgICBidWZbcCsrXSA9IHRyYW5zcGFyZW50X2luZGV4OyAgLy8gVHJhbnNwYXJlbnQgY29sb3IgaW5kZXguXG4gICAgICBidWZbcCsrXSA9IDA7ICAvLyBCbG9jayBUZXJtaW5hdG9yLlxuICAgIH1cblxuICAgIC8vIC0gSW1hZ2UgRGVzY3JpcHRvclxuICAgIGJ1ZltwKytdID0gMHgyYzsgIC8vIEltYWdlIFNlcGVyYXRvci5cbiAgICBidWZbcCsrXSA9IHggJiAweGZmOyBidWZbcCsrXSA9IHggPj4gOCAmIDB4ZmY7ICAvLyBMZWZ0LlxuICAgIGJ1ZltwKytdID0geSAmIDB4ZmY7IGJ1ZltwKytdID0geSA+PiA4ICYgMHhmZjsgIC8vIFRvcC5cbiAgICBidWZbcCsrXSA9IHcgJiAweGZmOyBidWZbcCsrXSA9IHcgPj4gOCAmIDB4ZmY7XG4gICAgYnVmW3ArK10gPSBoICYgMHhmZjsgYnVmW3ArK10gPSBoID4+IDggJiAweGZmO1xuICAgIC8vIE5PVEU6IE5vIHNvcnQgZmxhZyAodW51c2VkPykuXG4gICAgLy8gVE9ETyhkZWFubSk6IFN1cHBvcnQgaW50ZXJsYWNlLlxuICAgIGJ1ZltwKytdID0gdXNpbmdfbG9jYWxfcGFsZXR0ZSA9PT0gdHJ1ZSA/ICgweDgwIHwgKG1pbl9jb2RlX3NpemUtMSkpIDogMDtcblxuICAgIC8vIC0gTG9jYWwgQ29sb3IgVGFibGVcbiAgICBpZiAodXNpbmdfbG9jYWxfcGFsZXR0ZSA9PT0gdHJ1ZSkge1xuICAgICAgZm9yICh2YXIgaSA9IDAsIGlsID0gcGFsZXR0ZS5sZW5ndGg7IGkgPCBpbDsgKytpKSB7XG4gICAgICAgIHZhciByZ2IgPSBwYWxldHRlW2ldO1xuICAgICAgICBidWZbcCsrXSA9IHJnYiA+PiAxNiAmIDB4ZmY7XG4gICAgICAgIGJ1ZltwKytdID0gcmdiID4+IDggJiAweGZmO1xuICAgICAgICBidWZbcCsrXSA9IHJnYiAmIDB4ZmY7XG4gICAgICB9XG4gICAgfVxuXG4gICAgcCA9IEdpZldyaXRlck91dHB1dExaV0NvZGVTdHJlYW0oXG4gICAgICAgICAgICBidWYsIHAsIG1pbl9jb2RlX3NpemUgPCAyID8gMiA6IG1pbl9jb2RlX3NpemUsIGluZGV4ZWRfcGl4ZWxzKTtcbiAgfTtcblxuICB0aGlzLmVuZCA9IGZ1bmN0aW9uKCkge1xuICAgIGlmIChlbmRlZCA9PT0gZmFsc2UpIHtcbiAgICAgIGJ1ZltwKytdID0gMHgzYjsgIC8vIFRyYWlsZXIuXG4gICAgICBlbmRlZCA9IHRydWU7XG4gICAgfVxuICAgIHJldHVybiBwO1xuICB9O1xufVxuXG4vLyBNYWluIGNvbXByZXNzaW9uIHJvdXRpbmUsIHBhbGV0dGUgaW5kZXhlcyAtPiBMWlcgY29kZSBzdHJlYW0uXG4vLyB8aW5kZXhfc3RyZWFtfCBtdXN0IGhhdmUgYXQgbGVhc3Qgb25lIGVudHJ5LlxuZnVuY3Rpb24gR2lmV3JpdGVyT3V0cHV0TFpXQ29kZVN0cmVhbShidWYsIHAsIG1pbl9jb2RlX3NpemUsIGluZGV4X3N0cmVhbSkge1xuICBidWZbcCsrXSA9IG1pbl9jb2RlX3NpemU7XG4gIHZhciBjdXJfc3ViYmxvY2sgPSBwKys7ICAvLyBQb2ludGluZyBhdCB0aGUgbGVuZ3RoIGZpZWxkLlxuXG4gIHZhciBjbGVhcl9jb2RlID0gMSA8PCBtaW5fY29kZV9zaXplO1xuICB2YXIgY29kZV9tYXNrID0gY2xlYXJfY29kZSAtIDE7XG4gIHZhciBlb2lfY29kZSA9IGNsZWFyX2NvZGUgKyAxO1xuICB2YXIgbmV4dF9jb2RlID0gZW9pX2NvZGUgKyAxO1xuXG4gIHZhciBjdXJfY29kZV9zaXplID0gbWluX2NvZGVfc2l6ZSArIDE7ICAvLyBOdW1iZXIgb2YgYml0cyBwZXIgY29kZS5cbiAgdmFyIGN1cl9zaGlmdCA9IDA7XG4gIC8vIFdlIGhhdmUgYXQgbW9zdCAxMi1iaXQgY29kZXMsIHNvIHdlIHNob3VsZCBoYXZlIHRvIGhvbGQgYSBtYXggb2YgMTlcbiAgLy8gYml0cyBoZXJlIChhbmQgdGhlbiB3ZSB3b3VsZCB3cml0ZSBvdXQpLlxuICB2YXIgY3VyID0gMDtcblxuICBmdW5jdGlvbiBlbWl0X2J5dGVzX3RvX2J1ZmZlcihiaXRfYmxvY2tfc2l6ZSkge1xuICAgIHdoaWxlIChjdXJfc2hpZnQgPj0gYml0X2Jsb2NrX3NpemUpIHtcbiAgICAgIGJ1ZltwKytdID0gY3VyICYgMHhmZjtcbiAgICAgIGN1ciA+Pj0gODsgY3VyX3NoaWZ0IC09IDg7XG4gICAgICBpZiAocCA9PT0gY3VyX3N1YmJsb2NrICsgMjU2KSB7ICAvLyBGaW5pc2hlZCBhIHN1YmJsb2NrLlxuICAgICAgICBidWZbY3VyX3N1YmJsb2NrXSA9IDI1NTtcbiAgICAgICAgY3VyX3N1YmJsb2NrID0gcCsrO1xuICAgICAgfVxuICAgIH1cbiAgfVxuXG4gIGZ1bmN0aW9uIGVtaXRfY29kZShjKSB7XG4gICAgY3VyIHw9IGMgPDwgY3VyX3NoaWZ0O1xuICAgIGN1cl9zaGlmdCArPSBjdXJfY29kZV9zaXplO1xuICAgIGVtaXRfYnl0ZXNfdG9fYnVmZmVyKDgpO1xuICB9XG5cbiAgLy8gSSBhbSBub3QgYW4gZXhwZXJ0IG9uIHRoZSB0b3BpYywgYW5kIEkgZG9uJ3Qgd2FudCB0byB3cml0ZSBhIHRoZXNpcy5cbiAgLy8gSG93ZXZlciwgaXQgaXMgZ29vZCB0byBvdXRsaW5lIGhlcmUgdGhlIGJhc2ljIGFsZ29yaXRobSBhbmQgdGhlIGZldyBkYXRhXG4gIC8vIHN0cnVjdHVyZXMgYW5kIG9wdGltaXphdGlvbnMgaGVyZSB0aGF0IG1ha2UgdGhpcyBpbXBsZW1lbnRhdGlvbiBmYXN0LlxuICAvLyBUaGUgYmFzaWMgaWRlYSBiZWhpbmQgTFpXIGlzIHRvIGJ1aWxkIGEgdGFibGUgb2YgcHJldmlvdXNseSBzZWVuIHJ1bnNcbiAgLy8gYWRkcmVzc2VkIGJ5IGEgc2hvcnQgaWQgKGhlcmVpbiBjYWxsZWQgb3V0cHV0IGNvZGUpLiAgQWxsIGRhdGEgaXNcbiAgLy8gcmVmZXJlbmNlZCBieSBhIGNvZGUsIHdoaWNoIHJlcHJlc2VudHMgb25lIG9yIG1vcmUgdmFsdWVzIGZyb20gdGhlXG4gIC8vIG9yaWdpbmFsIGlucHV0IHN0cmVhbS4gIEFsbCBpbnB1dCBieXRlcyBjYW4gYmUgcmVmZXJlbmNlZCBhcyB0aGUgc2FtZVxuICAvLyB2YWx1ZSBhcyBhbiBvdXRwdXQgY29kZS4gIFNvIGlmIHlvdSBkaWRuJ3Qgd2FudCBhbnkgY29tcHJlc3Npb24sIHlvdVxuICAvLyBjb3VsZCBtb3JlIG9yIGxlc3MganVzdCBvdXRwdXQgdGhlIG9yaWdpbmFsIGJ5dGVzIGFzIGNvZGVzICh0aGVyZSBhcmVcbiAgLy8gc29tZSBkZXRhaWxzIHRvIHRoaXMsIGJ1dCBpdCBpcyB0aGUgaWRlYSkuICBJbiBvcmRlciB0byBhY2hpZXZlXG4gIC8vIGNvbXByZXNzaW9uLCB2YWx1ZXMgZ3JlYXRlciB0aGVuIHRoZSBpbnB1dCByYW5nZSAoY29kZXMgY2FuIGJlIHVwIHRvXG4gIC8vIDEyLWJpdCB3aGlsZSBpbnB1dCBvbmx5IDgtYml0KSByZXByZXNlbnQgYSBzZXF1ZW5jZSBvZiBwcmV2aW91c2x5IHNlZW5cbiAgLy8gaW5wdXRzLiAgVGhlIGRlY29tcHJlc3NvciBpcyBhYmxlIHRvIGJ1aWxkIHRoZSBzYW1lIG1hcHBpbmcgd2hpbGVcbiAgLy8gZGVjb2RpbmcsIHNvIHRoZXJlIGlzIGFsd2F5cyBhIHNoYXJlZCBjb21tb24ga25vd2xlZGdlIGJldHdlZW4gdGhlXG4gIC8vIGVuY29kaW5nIGFuZCBkZWNvZGVyLCB3aGljaCBpcyBhbHNvIGltcG9ydGFudCBmb3IgXCJ0aW1pbmdcIiBhc3BlY3RzIGxpa2VcbiAgLy8gaG93IHRvIGhhbmRsZSB2YXJpYWJsZSBiaXQgd2lkdGggY29kZSBlbmNvZGluZy5cbiAgLy9cbiAgLy8gT25lIG9idmlvdXMgYnV0IHZlcnkgaW1wb3J0YW50IGNvbnNlcXVlbmNlIG9mIHRoZSB0YWJsZSBzeXN0ZW0gaXMgdGhlcmVcbiAgLy8gaXMgYWx3YXlzIGEgdW5pcXVlIGlkIChhdCBtb3N0IDEyLWJpdHMpIHRvIG1hcCB0aGUgcnVucy4gICdBJyBtaWdodCBiZVxuICAvLyA0LCB0aGVuICdBQScgbWlnaHQgYmUgMTAsICdBQUEnIDExLCAnQUFBQScgMTIsIGV0Yy4gIFRoaXMgcmVsYXRpb25zaGlwXG4gIC8vIGNhbiBiZSB1c2VkIGZvciBhbiBlZmZlY2llbnQgbG9va3VwIHN0cmF0ZWd5IGZvciB0aGUgY29kZSBtYXBwaW5nLiAgV2VcbiAgLy8gbmVlZCB0byBrbm93IGlmIGEgcnVuIGhhcyBiZWVuIHNlZW4gYmVmb3JlLCBhbmQgYmUgYWJsZSB0byBtYXAgdGhhdCBydW5cbiAgLy8gdG8gdGhlIG91dHB1dCBjb2RlLiAgU2luY2Ugd2Ugc3RhcnQgd2l0aCBrbm93biB1bmlxdWUgaWRzIChpbnB1dCBieXRlcyksXG4gIC8vIGFuZCB0aGVuIGZyb20gdGhvc2UgYnVpbGQgbW9yZSB1bmlxdWUgaWRzICh0YWJsZSBlbnRyaWVzKSwgd2UgY2FuXG4gIC8vIGNvbnRpbnVlIHRoaXMgY2hhaW4gKGFsbW9zdCBsaWtlIGEgbGlua2VkIGxpc3QpIHRvIGFsd2F5cyBoYXZlIHNtYWxsXG4gIC8vIGludGVnZXIgdmFsdWVzIHRoYXQgcmVwcmVzZW50IHRoZSBjdXJyZW50IGJ5dGUgY2hhaW5zIGluIHRoZSBlbmNvZGVyLlxuICAvLyBUaGlzIG1lYW5zIGluc3RlYWQgb2YgdHJhY2tpbmcgdGhlIGlucHV0IGJ5dGVzIChBQUFBQkNEKSB0byBrbm93IG91clxuICAvLyBjdXJyZW50IHN0YXRlLCB3ZSBjYW4gdHJhY2sgdGhlIHRhYmxlIGVudHJ5IGZvciBBQUFBQkMgKGl0IGlzIGd1YXJhbnRlZWRcbiAgLy8gdG8gZXhpc3QgYnkgdGhlIG5hdHVyZSBvZiB0aGUgYWxnb3JpdGhtKSBhbmQgdGhlIG5leHQgY2hhcmFjdGVyIEQuXG4gIC8vIFRoZXJlZm9yIHRoZSB0dXBsZSBvZiAodGFibGVfZW50cnksIGJ5dGUpIGlzIGd1YXJhbnRlZWQgdG8gYWxzbyBiZVxuICAvLyB1bmlxdWUuICBUaGlzIGFsbG93cyB1cyB0byBjcmVhdGUgYSBzaW1wbGUgbG9va3VwIGtleSBmb3IgbWFwcGluZyBpbnB1dFxuICAvLyBzZXF1ZW5jZXMgdG8gY29kZXMgKHRhYmxlIGluZGljZXMpIHdpdGhvdXQgaGF2aW5nIHRvIHN0b3JlIG9yIHNlYXJjaFxuICAvLyBhbnkgb2YgdGhlIGNvZGUgc2VxdWVuY2VzLiAgU28gaWYgJ0FBQUEnIGhhcyBhIHRhYmxlIGVudHJ5IG9mIDEyLCB0aGVcbiAgLy8gdHVwbGUgb2YgKCdBQUFBJywgSykgZm9yIGFueSBpbnB1dCBieXRlIEsgd2lsbCBiZSB1bmlxdWUsIGFuZCBjYW4gYmUgb3VyXG4gIC8vIGtleS4gIFRoaXMgbGVhZHMgdG8gYSBpbnRlZ2VyIHZhbHVlIGF0IG1vc3QgMjAtYml0cywgd2hpY2ggY2FuIGFsd2F5c1xuICAvLyBmaXQgaW4gYW4gU01JIHZhbHVlIGFuZCBiZSB1c2VkIGFzIGEgZmFzdCBzcGFyc2UgYXJyYXkgLyBvYmplY3Qga2V5LlxuXG4gIC8vIE91dHB1dCBjb2RlIGZvciB0aGUgY3VycmVudCBjb250ZW50cyBvZiB0aGUgaW5kZXggYnVmZmVyLlxuICB2YXIgaWJfY29kZSA9IGluZGV4X3N0cmVhbVswXSAmIGNvZGVfbWFzazsgIC8vIExvYWQgZmlyc3QgaW5wdXQgaW5kZXguXG4gIHZhciBjb2RlX3RhYmxlID0geyB9OyAgLy8gS2V5J2Qgb24gb3VyIDIwLWJpdCBcInR1cGxlXCIuXG5cbiAgZW1pdF9jb2RlKGNsZWFyX2NvZGUpOyAgLy8gU3BlYyBzYXlzIGZpcnN0IGNvZGUgc2hvdWxkIGJlIGEgY2xlYXIgY29kZS5cblxuICAvLyBGaXJzdCBpbmRleCBhbHJlYWR5IGxvYWRlZCwgcHJvY2VzcyB0aGUgcmVzdCBvZiB0aGUgc3RyZWFtLlxuICBmb3IgKHZhciBpID0gMSwgaWwgPSBpbmRleF9zdHJlYW0ubGVuZ3RoOyBpIDwgaWw7ICsraSkge1xuICAgIHZhciBrID0gaW5kZXhfc3RyZWFtW2ldICYgY29kZV9tYXNrO1xuICAgIHZhciBjdXJfa2V5ID0gaWJfY29kZSA8PCA4IHwgazsgIC8vIChwcmV2LCBrKSB1bmlxdWUgdHVwbGUuXG4gICAgdmFyIGN1cl9jb2RlID0gY29kZV90YWJsZVtjdXJfa2V5XTsgIC8vIGJ1ZmZlciArIGsuXG5cbiAgICAvLyBDaGVjayBpZiB3ZSBoYXZlIHRvIGNyZWF0ZSBhIG5ldyBjb2RlIHRhYmxlIGVudHJ5LlxuICAgIGlmIChjdXJfY29kZSA9PT0gdW5kZWZpbmVkKSB7ICAvLyBXZSBkb24ndCBoYXZlIGJ1ZmZlciArIGsuXG4gICAgICAvLyBFbWl0IGluZGV4IGJ1ZmZlciAod2l0aG91dCBrKS5cbiAgICAgIC8vIFRoaXMgaXMgYW4gaW5saW5lIHZlcnNpb24gb2YgZW1pdF9jb2RlLCBiZWNhdXNlIHRoaXMgaXMgdGhlIGNvcmVcbiAgICAgIC8vIHdyaXRpbmcgcm91dGluZSBvZiB0aGUgY29tcHJlc3NvciAoYW5kIFY4IGNhbm5vdCBpbmxpbmUgZW1pdF9jb2RlXG4gICAgICAvLyBiZWNhdXNlIGl0IGlzIGEgY2xvc3VyZSBoZXJlIGluIGEgZGlmZmVyZW50IGNvbnRleHQpLiAgQWRkaXRpb25hbGx5XG4gICAgICAvLyB3ZSBjYW4gY2FsbCBlbWl0X2J5dGVfdG9fYnVmZmVyIGxlc3Mgb2Z0ZW4sIGJlY2F1c2Ugd2UgY2FuIGhhdmVcbiAgICAgIC8vIDMwLWJpdHMgKGZyb20gb3VyIDMxLWJpdCBzaWduZWQgU01JKSwgYW5kIHdlIGtub3cgb3VyIGNvZGVzIHdpbGwgb25seVxuICAgICAgLy8gYmUgMTItYml0cywgc28gY2FuIHNhZmVseSBoYXZlIDE4LWJpdHMgdGhlcmUgd2l0aG91dCBvdmVyZmxvdy5cbiAgICAgIC8vIGVtaXRfY29kZShpYl9jb2RlKTtcbiAgICAgIGN1ciB8PSBpYl9jb2RlIDw8IGN1cl9zaGlmdDtcbiAgICAgIGN1cl9zaGlmdCArPSBjdXJfY29kZV9zaXplO1xuICAgICAgd2hpbGUgKGN1cl9zaGlmdCA+PSA4KSB7XG4gICAgICAgIGJ1ZltwKytdID0gY3VyICYgMHhmZjtcbiAgICAgICAgY3VyID4+PSA4OyBjdXJfc2hpZnQgLT0gODtcbiAgICAgICAgaWYgKHAgPT09IGN1cl9zdWJibG9jayArIDI1NikgeyAgLy8gRmluaXNoZWQgYSBzdWJibG9jay5cbiAgICAgICAgICBidWZbY3VyX3N1YmJsb2NrXSA9IDI1NTtcbiAgICAgICAgICBjdXJfc3ViYmxvY2sgPSBwKys7XG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgaWYgKG5leHRfY29kZSA9PT0gNDA5NikgeyAgLy8gVGFibGUgZnVsbCwgbmVlZCBhIGNsZWFyLlxuICAgICAgICBlbWl0X2NvZGUoY2xlYXJfY29kZSk7XG4gICAgICAgIG5leHRfY29kZSA9IGVvaV9jb2RlICsgMTtcbiAgICAgICAgY3VyX2NvZGVfc2l6ZSA9IG1pbl9jb2RlX3NpemUgKyAxO1xuICAgICAgICBjb2RlX3RhYmxlID0geyB9O1xuICAgICAgfSBlbHNlIHsgIC8vIFRhYmxlIG5vdCBmdWxsLCBpbnNlcnQgYSBuZXcgZW50cnkuXG4gICAgICAgIC8vIEluY3JlYXNlIG91ciB2YXJpYWJsZSBiaXQgY29kZSBzaXplcyBpZiBuZWNlc3NhcnkuICBUaGlzIGlzIGEgYml0XG4gICAgICAgIC8vIHRyaWNreSBhcyBpdCBpcyBiYXNlZCBvbiBcInRpbWluZ1wiIGJldHdlZW4gdGhlIGVuY29kaW5nIGFuZFxuICAgICAgICAvLyBkZWNvZGVyLiAgRnJvbSB0aGUgZW5jb2RlcnMgcGVyc3BlY3RpdmUgdGhpcyBzaG91bGQgaGFwcGVuIGFmdGVyXG4gICAgICAgIC8vIHdlJ3ZlIGFscmVhZHkgZW1pdHRlZCB0aGUgaW5kZXggYnVmZmVyIGFuZCBhcmUgYWJvdXQgdG8gY3JlYXRlIHRoZVxuICAgICAgICAvLyBmaXJzdCB0YWJsZSBlbnRyeSB0aGF0IHdvdWxkIG92ZXJmbG93IG91ciBjdXJyZW50IGNvZGUgYml0IHNpemUuXG4gICAgICAgIGlmIChuZXh0X2NvZGUgPj0gKDEgPDwgY3VyX2NvZGVfc2l6ZSkpICsrY3VyX2NvZGVfc2l6ZTtcbiAgICAgICAgY29kZV90YWJsZVtjdXJfa2V5XSA9IG5leHRfY29kZSsrOyAgLy8gSW5zZXJ0IGludG8gY29kZSB0YWJsZS5cbiAgICAgIH1cblxuICAgICAgaWJfY29kZSA9IGs7ICAvLyBJbmRleCBidWZmZXIgdG8gc2luZ2xlIGlucHV0IGsuXG4gICAgfSBlbHNlIHtcbiAgICAgIGliX2NvZGUgPSBjdXJfY29kZTsgIC8vIEluZGV4IGJ1ZmZlciB0byBzZXF1ZW5jZSBpbiBjb2RlIHRhYmxlLlxuICAgIH1cbiAgfVxuXG4gIGVtaXRfY29kZShpYl9jb2RlKTsgIC8vIFRoZXJlIHdpbGwgc3RpbGwgYmUgc29tZXRoaW5nIGluIHRoZSBpbmRleCBidWZmZXIuXG4gIGVtaXRfY29kZShlb2lfY29kZSk7ICAvLyBFbmQgT2YgSW5mb3JtYXRpb24uXG5cbiAgLy8gRmx1c2ggLyBmaW5hbGl6ZSB0aGUgc3ViLWJsb2NrcyBzdHJlYW0gdG8gdGhlIGJ1ZmZlci5cbiAgZW1pdF9ieXRlc190b19idWZmZXIoMSk7XG5cbiAgLy8gRmluaXNoIHRoZSBzdWItYmxvY2tzLCB3cml0aW5nIG91dCBhbnkgdW5maW5pc2hlZCBsZW5ndGhzIGFuZFxuICAvLyB0ZXJtaW5hdGluZyB3aXRoIGEgc3ViLWJsb2NrIG9mIGxlbmd0aCAwLiAgSWYgd2UgaGF2ZSBhbHJlYWR5IHN0YXJ0ZWRcbiAgLy8gYnV0IG5vdCB5ZXQgdXNlZCBhIHN1Yi1ibG9jayBpdCBjYW4ganVzdCBiZWNvbWUgdGhlIHRlcm1pbmF0b3IuXG4gIGlmIChjdXJfc3ViYmxvY2sgKyAxID09PSBwKSB7ICAvLyBTdGFydGVkIGJ1dCB1bnVzZWQuXG4gICAgYnVmW2N1cl9zdWJibG9ja10gPSAwO1xuICB9IGVsc2UgeyAgLy8gU3RhcnRlZCBhbmQgdXNlZCwgd3JpdGUgbGVuZ3RoIGFuZCBhZGRpdGlvbmFsIHRlcm1pbmF0b3IgYmxvY2suXG4gICAgYnVmW2N1cl9zdWJibG9ja10gPSBwIC0gY3VyX3N1YmJsb2NrIC0gMTtcbiAgICBidWZbcCsrXSA9IDA7XG4gIH1cbiAgcmV0dXJuIHA7XG59XG5cbmZ1bmN0aW9uIEdpZlJlYWRlcihidWYpIHtcbiAgdmFyIHAgPSAwO1xuXG4gIC8vIC0gSGVhZGVyIChHSUY4N2Egb3IgR0lGODlhKS5cbiAgaWYgKGJ1ZltwKytdICE9PSAweDQ3IHx8ICAgICAgICAgICAgYnVmW3ArK10gIT09IDB4NDkgfHwgYnVmW3ArK10gIT09IDB4NDYgfHxcbiAgICAgIGJ1ZltwKytdICE9PSAweDM4IHx8IChidWZbcCsrXSsxICYgMHhmZCkgIT09IDB4MzggfHwgYnVmW3ArK10gIT09IDB4NjEpIHtcbiAgICB0aHJvdyBcIkludmFsaWQgR0lGIDg3YS84OWEgaGVhZGVyLlwiO1xuICB9XG5cbiAgLy8gLSBMb2dpY2FsIFNjcmVlbiBEZXNjcmlwdG9yLlxuICB2YXIgd2lkdGggPSBidWZbcCsrXSB8IGJ1ZltwKytdIDw8IDg7XG4gIHZhciBoZWlnaHQgPSBidWZbcCsrXSB8IGJ1ZltwKytdIDw8IDg7XG4gIHZhciBwZjAgPSBidWZbcCsrXTsgIC8vIDxQYWNrZWQgRmllbGRzPi5cbiAgdmFyIGdsb2JhbF9wYWxldHRlX2ZsYWcgPSBwZjAgPj4gNztcbiAgdmFyIG51bV9nbG9iYWxfY29sb3JzX3BvdzIgPSBwZjAgJiAweDc7XG4gIHZhciBudW1fZ2xvYmFsX2NvbG9ycyA9IDEgPDwgKG51bV9nbG9iYWxfY29sb3JzX3BvdzIgKyAxKTtcbiAgdmFyIGJhY2tncm91bmQgPSBidWZbcCsrXTtcbiAgYnVmW3ArK107ICAvLyBQaXhlbCBhc3BlY3QgcmF0aW8gKHVudXNlZD8pLlxuXG4gIHZhciBnbG9iYWxfcGFsZXR0ZV9vZmZzZXQgPSBudWxsO1xuXG4gIGlmIChnbG9iYWxfcGFsZXR0ZV9mbGFnKSB7XG4gICAgZ2xvYmFsX3BhbGV0dGVfb2Zmc2V0ID0gcDtcbiAgICBwICs9IG51bV9nbG9iYWxfY29sb3JzICogMzsgIC8vIFNlZWsgcGFzdCBwYWxldHRlLlxuICB9XG5cbiAgdmFyIG5vX2VvZiA9IHRydWU7XG5cbiAgdmFyIGZyYW1lcyA9IFsgXTtcblxuICB2YXIgZGVsYXkgPSAwO1xuICB2YXIgdHJhbnNwYXJlbnRfaW5kZXggPSBudWxsO1xuICB2YXIgZGlzcG9zYWwgPSAwOyAgLy8gMCAtIE5vIGRpc3Bvc2FsIHNwZWNpZmllZC5cbiAgdmFyIGxvb3BfY291bnQgPSBudWxsO1xuXG4gIHRoaXMud2lkdGggPSB3aWR0aDtcbiAgdGhpcy5oZWlnaHQgPSBoZWlnaHQ7XG5cbiAgd2hpbGUgKG5vX2VvZiAmJiBwIDwgYnVmLmxlbmd0aCkge1xuICAgIHN3aXRjaCAoYnVmW3ArK10pIHtcbiAgICAgIGNhc2UgMHgyMTogIC8vIEdyYXBoaWNzIENvbnRyb2wgRXh0ZW5zaW9uIEJsb2NrXG4gICAgICAgIHN3aXRjaCAoYnVmW3ArK10pIHtcbiAgICAgICAgICBjYXNlIDB4ZmY6ICAvLyBBcHBsaWNhdGlvbiBzcGVjaWZpYyBibG9ja1xuICAgICAgICAgICAgLy8gVHJ5IGlmIGl0J3MgYSBOZXRzY2FwZSBibG9jayAod2l0aCBhbmltYXRpb24gbG9vcCBjb3VudGVyKS5cbiAgICAgICAgICAgIGlmIChidWZbcCAgIF0gIT09IDB4MGIgfHwgIC8vIDIxIEZGIGFscmVhZHkgcmVhZCwgY2hlY2sgYmxvY2sgc2l6ZS5cbiAgICAgICAgICAgICAgICAvLyBORVRTQ0FQRTIuMFxuICAgICAgICAgICAgICAgIGJ1ZltwKzEgXSA9PSAweDRlICYmIGJ1ZltwKzIgXSA9PSAweDQ1ICYmIGJ1ZltwKzMgXSA9PSAweDU0ICYmXG4gICAgICAgICAgICAgICAgYnVmW3ArNCBdID09IDB4NTMgJiYgYnVmW3ArNSBdID09IDB4NDMgJiYgYnVmW3ArNiBdID09IDB4NDEgJiZcbiAgICAgICAgICAgICAgICBidWZbcCs3IF0gPT0gMHg1MCAmJiBidWZbcCs4IF0gPT0gMHg0NSAmJiBidWZbcCs5IF0gPT0gMHgzMiAmJlxuICAgICAgICAgICAgICAgIGJ1ZltwKzEwXSA9PSAweDJlICYmIGJ1ZltwKzExXSA9PSAweDMwICYmXG4gICAgICAgICAgICAgICAgLy8gU3ViLWJsb2NrXG4gICAgICAgICAgICAgICAgYnVmW3ArMTJdID09IDB4MDMgJiYgYnVmW3ArMTNdID09IDB4MDEgJiYgYnVmW3ArMTZdID09IDApIHtcbiAgICAgICAgICAgICAgcCArPSAxNDtcbiAgICAgICAgICAgICAgbG9vcF9jb3VudCA9IGJ1ZltwKytdIHwgYnVmW3ArK10gPDwgODtcbiAgICAgICAgICAgICAgcCsrOyAgLy8gU2tpcCB0ZXJtaW5hdG9yLlxuICAgICAgICAgICAgfSBlbHNlIHsgIC8vIFdlIGRvbid0IGtub3cgd2hhdCBpdCBpcywganVzdCB0cnkgdG8gZ2V0IHBhc3QgaXQuXG4gICAgICAgICAgICAgIHAgKz0gMTI7XG4gICAgICAgICAgICAgIHdoaWxlICh0cnVlKSB7ICAvLyBTZWVrIHRocm91Z2ggc3ViYmxvY2tzLlxuICAgICAgICAgICAgICAgIHZhciBibG9ja19zaXplID0gYnVmW3ArK107XG4gICAgICAgICAgICAgICAgaWYgKGJsb2NrX3NpemUgPT09IDApIGJyZWFrO1xuICAgICAgICAgICAgICAgIHAgKz0gYmxvY2tfc2l6ZTtcbiAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICAgYnJlYWs7XG5cbiAgICAgICAgICBjYXNlIDB4Zjk6ICAvLyBHcmFwaGljcyBDb250cm9sIEV4dGVuc2lvblxuICAgICAgICAgICAgaWYgKGJ1ZltwKytdICE9PSAweDQgfHwgYnVmW3ArNF0gIT09IDApXG4gICAgICAgICAgICAgIHRocm93IFwiSW52YWxpZCBncmFwaGljcyBleHRlbnNpb24gYmxvY2suXCI7XG4gICAgICAgICAgICB2YXIgcGYxID0gYnVmW3ArK107XG4gICAgICAgICAgICBkZWxheSA9IGJ1ZltwKytdIHwgYnVmW3ArK10gPDwgODtcbiAgICAgICAgICAgIHRyYW5zcGFyZW50X2luZGV4ID0gYnVmW3ArK107XG4gICAgICAgICAgICBpZiAoKHBmMSAmIDEpID09PSAwKSB0cmFuc3BhcmVudF9pbmRleCA9IG51bGw7XG4gICAgICAgICAgICBkaXNwb3NhbCA9IHBmMSA+PiAyICYgMHg3O1xuICAgICAgICAgICAgcCsrOyAgLy8gU2tpcCB0ZXJtaW5hdG9yLlxuICAgICAgICAgICAgYnJlYWs7XG5cbiAgICAgICAgICBjYXNlIDB4ZmU6ICAvLyBDb21tZW50IEV4dGVuc2lvbi5cbiAgICAgICAgICAgIHdoaWxlICh0cnVlKSB7ICAvLyBTZWVrIHRocm91Z2ggc3ViYmxvY2tzLlxuICAgICAgICAgICAgICB2YXIgYmxvY2tfc2l6ZSA9IGJ1ZltwKytdO1xuICAgICAgICAgICAgICBpZiAoYmxvY2tfc2l6ZSA9PT0gMCkgYnJlYWs7XG4gICAgICAgICAgICAgIC8vIGNvbnNvbGUubG9nKGJ1Zi5zbGljZShwLCBwK2Jsb2NrX3NpemUpLnRvU3RyaW5nKCdhc2NpaScpKTtcbiAgICAgICAgICAgICAgcCArPSBibG9ja19zaXplO1xuICAgICAgICAgICAgfVxuICAgICAgICAgICAgYnJlYWs7XG5cbiAgICAgICAgICBkZWZhdWx0OlxuICAgICAgICAgICAgdGhyb3cgXCJVbmtub3duIGdyYXBoaWMgY29udHJvbCBsYWJlbDogMHhcIiArIGJ1ZltwLTFdLnRvU3RyaW5nKDE2KTtcbiAgICAgICAgfVxuICAgICAgICBicmVhaztcblxuICAgICAgY2FzZSAweDJjOiAgLy8gSW1hZ2UgRGVzY3JpcHRvci5cbiAgICAgICAgdmFyIHggPSBidWZbcCsrXSB8IGJ1ZltwKytdIDw8IDg7XG4gICAgICAgIHZhciB5ID0gYnVmW3ArK10gfCBidWZbcCsrXSA8PCA4O1xuICAgICAgICB2YXIgdyA9IGJ1ZltwKytdIHwgYnVmW3ArK10gPDwgODtcbiAgICAgICAgdmFyIGggPSBidWZbcCsrXSB8IGJ1ZltwKytdIDw8IDg7XG4gICAgICAgIHZhciBwZjIgPSBidWZbcCsrXTtcbiAgICAgICAgdmFyIGxvY2FsX3BhbGV0dGVfZmxhZyA9IHBmMiA+PiA3O1xuICAgICAgICB2YXIgaW50ZXJsYWNlX2ZsYWcgPSBwZjIgPj4gNiAmIDE7XG4gICAgICAgIHZhciBudW1fbG9jYWxfY29sb3JzX3BvdzIgPSBwZjIgJiAweDc7XG4gICAgICAgIHZhciBudW1fbG9jYWxfY29sb3JzID0gMSA8PCAobnVtX2xvY2FsX2NvbG9yc19wb3cyICsgMSk7XG4gICAgICAgIHZhciBwYWxldHRlX29mZnNldCA9IGdsb2JhbF9wYWxldHRlX29mZnNldDtcbiAgICAgICAgdmFyIGhhc19sb2NhbF9wYWxldHRlID0gZmFsc2U7XG4gICAgICAgIGlmIChsb2NhbF9wYWxldHRlX2ZsYWcpIHtcbiAgICAgICAgICB2YXIgaGFzX2xvY2FsX3BhbGV0dGUgPSB0cnVlO1xuICAgICAgICAgIHBhbGV0dGVfb2Zmc2V0ID0gcDsgIC8vIE92ZXJyaWRlIHdpdGggbG9jYWwgcGFsZXR0ZS5cbiAgICAgICAgICBwICs9IG51bV9sb2NhbF9jb2xvcnMgKiAzOyAgLy8gU2VlayBwYXN0IHBhbGV0dGUuXG4gICAgICAgIH1cblxuICAgICAgICB2YXIgZGF0YV9vZmZzZXQgPSBwO1xuXG4gICAgICAgIHArKzsgIC8vIGNvZGVzaXplXG4gICAgICAgIHdoaWxlICh0cnVlKSB7XG4gICAgICAgICAgdmFyIGJsb2NrX3NpemUgPSBidWZbcCsrXTtcbiAgICAgICAgICBpZiAoYmxvY2tfc2l6ZSA9PT0gMCkgYnJlYWs7XG4gICAgICAgICAgcCArPSBibG9ja19zaXplO1xuICAgICAgICB9XG5cbiAgICAgICAgZnJhbWVzLnB1c2goe3g6IHgsIHk6IHksIHdpZHRoOiB3LCBoZWlnaHQ6IGgsXG4gICAgICAgICAgICAgICAgICAgICBoYXNfbG9jYWxfcGFsZXR0ZTogaGFzX2xvY2FsX3BhbGV0dGUsXG4gICAgICAgICAgICAgICAgICAgICBwYWxldHRlX29mZnNldDogcGFsZXR0ZV9vZmZzZXQsXG4gICAgICAgICAgICAgICAgICAgICBkYXRhX29mZnNldDogZGF0YV9vZmZzZXQsXG4gICAgICAgICAgICAgICAgICAgICBkYXRhX2xlbmd0aDogcCAtIGRhdGFfb2Zmc2V0LFxuICAgICAgICAgICAgICAgICAgICAgdHJhbnNwYXJlbnRfaW5kZXg6IHRyYW5zcGFyZW50X2luZGV4LFxuICAgICAgICAgICAgICAgICAgICAgaW50ZXJsYWNlZDogISFpbnRlcmxhY2VfZmxhZyxcbiAgICAgICAgICAgICAgICAgICAgIGRlbGF5OiBkZWxheSxcbiAgICAgICAgICAgICAgICAgICAgIGRpc3Bvc2FsOiBkaXNwb3NhbH0pO1xuICAgICAgICBicmVhaztcblxuICAgICAgY2FzZSAweDNiOiAgLy8gVHJhaWxlciBNYXJrZXIgKGVuZCBvZiBmaWxlKS5cbiAgICAgICAgbm9fZW9mID0gZmFsc2U7XG4gICAgICAgIGJyZWFrO1xuXG4gICAgICBkZWZhdWx0OlxuICAgICAgICB0aHJvdyBcIlVua25vd24gZ2lmIGJsb2NrOiAweFwiICsgYnVmW3AtMV0udG9TdHJpbmcoMTYpO1xuICAgICAgICBicmVhaztcbiAgICB9XG4gIH1cblxuICB0aGlzLm51bUZyYW1lcyA9IGZ1bmN0aW9uKCkge1xuICAgIHJldHVybiBmcmFtZXMubGVuZ3RoO1xuICB9O1xuXG4gIHRoaXMubG9vcENvdW50ID0gZnVuY3Rpb24oKSB7XG4gICAgcmV0dXJuIGxvb3BfY291bnQ7XG4gIH07XG5cbiAgdGhpcy5mcmFtZUluZm8gPSBmdW5jdGlvbihmcmFtZV9udW0pIHtcbiAgICBpZiAoZnJhbWVfbnVtIDwgMCB8fCBmcmFtZV9udW0gPj0gZnJhbWVzLmxlbmd0aClcbiAgICAgIHRocm93IFwiRnJhbWUgaW5kZXggb3V0IG9mIHJhbmdlLlwiO1xuICAgIHJldHVybiBmcmFtZXNbZnJhbWVfbnVtXTtcbiAgfVxuXG4gIHRoaXMuZGVjb2RlQW5kQmxpdEZyYW1lQkdSQSA9IGZ1bmN0aW9uKGZyYW1lX251bSwgcGl4ZWxzKSB7XG4gICAgdmFyIGZyYW1lID0gdGhpcy5mcmFtZUluZm8oZnJhbWVfbnVtKTtcbiAgICB2YXIgbnVtX3BpeGVscyA9IGZyYW1lLndpZHRoICogZnJhbWUuaGVpZ2h0O1xuICAgIHZhciBpbmRleF9zdHJlYW0gPSBuZXcgVWludDhBcnJheShudW1fcGl4ZWxzKTsgIC8vIEF0IG1vc3QgOC1iaXQgaW5kaWNlcy5cbiAgICBHaWZSZWFkZXJMWldPdXRwdXRJbmRleFN0cmVhbShcbiAgICAgICAgYnVmLCBmcmFtZS5kYXRhX29mZnNldCwgaW5kZXhfc3RyZWFtLCBudW1fcGl4ZWxzKTtcbiAgICB2YXIgcGFsZXR0ZV9vZmZzZXQgPSBmcmFtZS5wYWxldHRlX29mZnNldDtcblxuICAgIC8vIE5PVEUoZGVhbm0pOiBJdCBzZWVtcyB0byBiZSBtdWNoIGZhc3RlciB0byBjb21wYXJlIGluZGV4IHRvIDI1NiB0aGFuXG4gICAgLy8gdG8gPT09IG51bGwuICBOb3Qgc3VyZSB3aHksIGJ1dCBDb21wYXJlU3R1Yl9FUV9TVFJJQ1Qgc2hvd3MgdXAgaGlnaCBpblxuICAgIC8vIHRoZSBwcm9maWxlLCBub3Qgc3VyZSBpZiBpdCdzIHJlbGF0ZWQgdG8gdXNpbmcgYSBVaW50OEFycmF5LlxuICAgIHZhciB0cmFucyA9IGZyYW1lLnRyYW5zcGFyZW50X2luZGV4O1xuICAgIGlmICh0cmFucyA9PT0gbnVsbCkgdHJhbnMgPSAyNTY7XG5cbiAgICAvLyBXZSBhcmUgcG9zc2libHkganVzdCBibGl0dGluZyB0byBhIHBvcnRpb24gb2YgdGhlIGVudGlyZSBmcmFtZS5cbiAgICAvLyBUaGF0IGlzIGEgc3VicmVjdCB3aXRoaW4gdGhlIGZyYW1lcmVjdCwgc28gdGhlIGFkZGl0aW9uYWwgcGl4ZWxzXG4gICAgLy8gbXVzdCBiZSBza2lwcGVkIG92ZXIgYWZ0ZXIgd2UgZmluaXNoZWQgYSBzY2FubGluZS5cbiAgICB2YXIgZnJhbWV3aWR0aCAgPSBmcmFtZS53aWR0aDtcbiAgICB2YXIgZnJhbWVzdHJpZGUgPSB3aWR0aCAtIGZyYW1ld2lkdGg7XG4gICAgdmFyIHhsZWZ0ICAgICAgID0gZnJhbWV3aWR0aDsgIC8vIE51bWJlciBvZiBzdWJyZWN0IHBpeGVscyBsZWZ0IGluIHNjYW5saW5lLlxuXG4gICAgLy8gT3V0cHV0IGluZGljaWVzIG9mIHRoZSB0b3AgbGVmdCBhbmQgYm90dG9tIHJpZ2h0IGNvcm5lcnMgb2YgdGhlIHN1YnJlY3QuXG4gICAgdmFyIG9wYmVnID0gKChmcmFtZS55ICogd2lkdGgpICsgZnJhbWUueCkgKiA0O1xuICAgIHZhciBvcGVuZCA9ICgoZnJhbWUueSArIGZyYW1lLmhlaWdodCkgKiB3aWR0aCArIGZyYW1lLngpICogNDtcbiAgICB2YXIgb3AgICAgPSBvcGJlZztcblxuICAgIHZhciBzY2Fuc3RyaWRlID0gZnJhbWVzdHJpZGUgKiA0O1xuXG4gICAgLy8gVXNlIHNjYW5zdHJpZGUgdG8gc2tpcCBwYXN0IHRoZSByb3dzIHdoZW4gaW50ZXJsYWNpbmcuICBUaGlzIGlzIHNraXBwaW5nXG4gICAgLy8gNyByb3dzIGZvciB0aGUgZmlyc3QgdHdvIHBhc3NlcywgdGhlbiAzIHRoZW4gMS5cbiAgICBpZiAoZnJhbWUuaW50ZXJsYWNlZCA9PT0gdHJ1ZSkge1xuICAgICAgc2NhbnN0cmlkZSArPSB3aWR0aCAqIDQgKiA3OyAgLy8gUGFzcyAxLlxuICAgIH1cblxuICAgIHZhciBpbnRlcmxhY2Vza2lwID0gODsgIC8vIFRyYWNraW5nIHRoZSByb3cgaW50ZXJ2YWwgaW4gdGhlIGN1cnJlbnQgcGFzcy5cblxuICAgIGZvciAodmFyIGkgPSAwLCBpbCA9IGluZGV4X3N0cmVhbS5sZW5ndGg7IGkgPCBpbDsgKytpKSB7XG4gICAgICB2YXIgaW5kZXggPSBpbmRleF9zdHJlYW1baV07XG5cbiAgICAgIGlmICh4bGVmdCA9PT0gMCkgeyAgLy8gQmVnaW5uaW5nIG9mIG5ldyBzY2FuIGxpbmVcbiAgICAgICAgb3AgKz0gc2NhbnN0cmlkZTtcbiAgICAgICAgeGxlZnQgPSBmcmFtZXdpZHRoO1xuICAgICAgICBpZiAob3AgPj0gb3BlbmQpIHsgLy8gQ2F0Y2ggdGhlIHdyYXAgdG8gc3dpdGNoIHBhc3NlcyB3aGVuIGludGVybGFjaW5nLlxuICAgICAgICAgIHNjYW5zdHJpZGUgPSBmcmFtZXN0cmlkZSAqIDQgKyB3aWR0aCAqIDQgKiAoaW50ZXJsYWNlc2tpcC0xKTtcbiAgICAgICAgICAvLyBpbnRlcmxhY2Vza2lwIC8gMiAqIDQgaXMgaW50ZXJsYWNlc2tpcCA8PCAxLlxuICAgICAgICAgIG9wID0gb3BiZWcgKyAoZnJhbWV3aWR0aCArIGZyYW1lc3RyaWRlKSAqIChpbnRlcmxhY2Vza2lwIDw8IDEpO1xuICAgICAgICAgIGludGVybGFjZXNraXAgPj49IDE7XG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgaWYgKGluZGV4ID09PSB0cmFucykge1xuICAgICAgICBvcCArPSA0O1xuICAgICAgfSBlbHNlIHtcbiAgICAgICAgdmFyIHIgPSBidWZbcGFsZXR0ZV9vZmZzZXQgKyBpbmRleCAqIDNdO1xuICAgICAgICB2YXIgZyA9IGJ1ZltwYWxldHRlX29mZnNldCArIGluZGV4ICogMyArIDFdO1xuICAgICAgICB2YXIgYiA9IGJ1ZltwYWxldHRlX29mZnNldCArIGluZGV4ICogMyArIDJdO1xuICAgICAgICBwaXhlbHNbb3ArK10gPSBiO1xuICAgICAgICBwaXhlbHNbb3ArK10gPSBnO1xuICAgICAgICBwaXhlbHNbb3ArK10gPSByO1xuICAgICAgICBwaXhlbHNbb3ArK10gPSAyNTU7XG4gICAgICB9XG4gICAgICAtLXhsZWZ0O1xuICAgIH1cbiAgfTtcblxuICAvLyBJIHdpbGwgZ28gdG8gY29weSBhbmQgcGFzdGUgaGVsbCBvbmUgZGF5Li4uXG4gIHRoaXMuZGVjb2RlQW5kQmxpdEZyYW1lUkdCQSA9IGZ1bmN0aW9uKGZyYW1lX251bSwgcGl4ZWxzKSB7XG4gICAgdmFyIGZyYW1lID0gdGhpcy5mcmFtZUluZm8oZnJhbWVfbnVtKTtcbiAgICB2YXIgbnVtX3BpeGVscyA9IGZyYW1lLndpZHRoICogZnJhbWUuaGVpZ2h0O1xuICAgIHZhciBpbmRleF9zdHJlYW0gPSBuZXcgVWludDhBcnJheShudW1fcGl4ZWxzKTsgIC8vIEF0IG1vc3QgOC1iaXQgaW5kaWNlcy5cbiAgICBHaWZSZWFkZXJMWldPdXRwdXRJbmRleFN0cmVhbShcbiAgICAgICAgYnVmLCBmcmFtZS5kYXRhX29mZnNldCwgaW5kZXhfc3RyZWFtLCBudW1fcGl4ZWxzKTtcbiAgICB2YXIgcGFsZXR0ZV9vZmZzZXQgPSBmcmFtZS5wYWxldHRlX29mZnNldDtcblxuICAgIC8vIE5PVEUoZGVhbm0pOiBJdCBzZWVtcyB0byBiZSBtdWNoIGZhc3RlciB0byBjb21wYXJlIGluZGV4IHRvIDI1NiB0aGFuXG4gICAgLy8gdG8gPT09IG51bGwuICBOb3Qgc3VyZSB3aHksIGJ1dCBDb21wYXJlU3R1Yl9FUV9TVFJJQ1Qgc2hvd3MgdXAgaGlnaCBpblxuICAgIC8vIHRoZSBwcm9maWxlLCBub3Qgc3VyZSBpZiBpdCdzIHJlbGF0ZWQgdG8gdXNpbmcgYSBVaW50OEFycmF5LlxuICAgIHZhciB0cmFucyA9IGZyYW1lLnRyYW5zcGFyZW50X2luZGV4O1xuICAgIGlmICh0cmFucyA9PT0gbnVsbCkgdHJhbnMgPSAyNTY7XG5cbiAgICAvLyBXZSBhcmUgcG9zc2libHkganVzdCBibGl0dGluZyB0byBhIHBvcnRpb24gb2YgdGhlIGVudGlyZSBmcmFtZS5cbiAgICAvLyBUaGF0IGlzIGEgc3VicmVjdCB3aXRoaW4gdGhlIGZyYW1lcmVjdCwgc28gdGhlIGFkZGl0aW9uYWwgcGl4ZWxzXG4gICAgLy8gbXVzdCBiZSBza2lwcGVkIG92ZXIgYWZ0ZXIgd2UgZmluaXNoZWQgYSBzY2FubGluZS5cbiAgICB2YXIgZnJhbWV3aWR0aCAgPSBmcmFtZS53aWR0aDtcbiAgICB2YXIgZnJhbWVzdHJpZGUgPSB3aWR0aCAtIGZyYW1ld2lkdGg7XG4gICAgdmFyIHhsZWZ0ICAgICAgID0gZnJhbWV3aWR0aDsgIC8vIE51bWJlciBvZiBzdWJyZWN0IHBpeGVscyBsZWZ0IGluIHNjYW5saW5lLlxuXG4gICAgLy8gT3V0cHV0IGluZGljaWVzIG9mIHRoZSB0b3AgbGVmdCBhbmQgYm90dG9tIHJpZ2h0IGNvcm5lcnMgb2YgdGhlIHN1YnJlY3QuXG4gICAgdmFyIG9wYmVnID0gKChmcmFtZS55ICogd2lkdGgpICsgZnJhbWUueCkgKiA0O1xuICAgIHZhciBvcGVuZCA9ICgoZnJhbWUueSArIGZyYW1lLmhlaWdodCkgKiB3aWR0aCArIGZyYW1lLngpICogNDtcbiAgICB2YXIgb3AgICAgPSBvcGJlZztcblxuICAgIHZhciBzY2Fuc3RyaWRlID0gZnJhbWVzdHJpZGUgKiA0O1xuXG4gICAgLy8gVXNlIHNjYW5zdHJpZGUgdG8gc2tpcCBwYXN0IHRoZSByb3dzIHdoZW4gaW50ZXJsYWNpbmcuICBUaGlzIGlzIHNraXBwaW5nXG4gICAgLy8gNyByb3dzIGZvciB0aGUgZmlyc3QgdHdvIHBhc3NlcywgdGhlbiAzIHRoZW4gMS5cbiAgICBpZiAoZnJhbWUuaW50ZXJsYWNlZCA9PT0gdHJ1ZSkge1xuICAgICAgc2NhbnN0cmlkZSArPSB3aWR0aCAqIDQgKiA3OyAgLy8gUGFzcyAxLlxuICAgIH1cblxuICAgIHZhciBpbnRlcmxhY2Vza2lwID0gODsgIC8vIFRyYWNraW5nIHRoZSByb3cgaW50ZXJ2YWwgaW4gdGhlIGN1cnJlbnQgcGFzcy5cblxuICAgIGZvciAodmFyIGkgPSAwLCBpbCA9IGluZGV4X3N0cmVhbS5sZW5ndGg7IGkgPCBpbDsgKytpKSB7XG4gICAgICB2YXIgaW5kZXggPSBpbmRleF9zdHJlYW1baV07XG5cbiAgICAgIGlmICh4bGVmdCA9PT0gMCkgeyAgLy8gQmVnaW5uaW5nIG9mIG5ldyBzY2FuIGxpbmVcbiAgICAgICAgb3AgKz0gc2NhbnN0cmlkZTtcbiAgICAgICAgeGxlZnQgPSBmcmFtZXdpZHRoO1xuICAgICAgICBpZiAob3AgPj0gb3BlbmQpIHsgLy8gQ2F0Y2ggdGhlIHdyYXAgdG8gc3dpdGNoIHBhc3NlcyB3aGVuIGludGVybGFjaW5nLlxuICAgICAgICAgIHNjYW5zdHJpZGUgPSBmcmFtZXN0cmlkZSAqIDQgKyB3aWR0aCAqIDQgKiAoaW50ZXJsYWNlc2tpcC0xKTtcbiAgICAgICAgICAvLyBpbnRlcmxhY2Vza2lwIC8gMiAqIDQgaXMgaW50ZXJsYWNlc2tpcCA8PCAxLlxuICAgICAgICAgIG9wID0gb3BiZWcgKyAoZnJhbWV3aWR0aCArIGZyYW1lc3RyaWRlKSAqIChpbnRlcmxhY2Vza2lwIDw8IDEpO1xuICAgICAgICAgIGludGVybGFjZXNraXAgPj49IDE7XG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgaWYgKGluZGV4ID09PSB0cmFucykge1xuICAgICAgICBvcCArPSA0O1xuICAgICAgfSBlbHNlIHtcbiAgICAgICAgdmFyIHIgPSBidWZbcGFsZXR0ZV9vZmZzZXQgKyBpbmRleCAqIDNdO1xuICAgICAgICB2YXIgZyA9IGJ1ZltwYWxldHRlX29mZnNldCArIGluZGV4ICogMyArIDFdO1xuICAgICAgICB2YXIgYiA9IGJ1ZltwYWxldHRlX29mZnNldCArIGluZGV4ICogMyArIDJdO1xuICAgICAgICBwaXhlbHNbb3ArK10gPSByO1xuICAgICAgICBwaXhlbHNbb3ArK10gPSBnO1xuICAgICAgICBwaXhlbHNbb3ArK10gPSBiO1xuICAgICAgICBwaXhlbHNbb3ArK10gPSAyNTU7XG4gICAgICB9XG4gICAgICAtLXhsZWZ0O1xuICAgIH1cbiAgfTtcbn1cblxuZnVuY3Rpb24gR2lmUmVhZGVyTFpXT3V0cHV0SW5kZXhTdHJlYW0oY29kZV9zdHJlYW0sIHAsIG91dHB1dCwgb3V0cHV0X2xlbmd0aCkge1xuICB2YXIgbWluX2NvZGVfc2l6ZSA9IGNvZGVfc3RyZWFtW3ArK107XG5cbiAgdmFyIGNsZWFyX2NvZGUgPSAxIDw8IG1pbl9jb2RlX3NpemU7XG4gIHZhciBlb2lfY29kZSA9IGNsZWFyX2NvZGUgKyAxO1xuICB2YXIgbmV4dF9jb2RlID0gZW9pX2NvZGUgKyAxO1xuXG4gIHZhciBjdXJfY29kZV9zaXplID0gbWluX2NvZGVfc2l6ZSArIDE7ICAvLyBOdW1iZXIgb2YgYml0cyBwZXIgY29kZS5cbiAgLy8gTk9URTogVGhpcyBzaGFyZXMgdGhlIHNhbWUgbmFtZSBhcyB0aGUgZW5jb2RlciwgYnV0IGhhcyBhIGRpZmZlcmVudFxuICAvLyBtZWFuaW5nIGhlcmUuICBIZXJlIHRoaXMgbWFza3MgZWFjaCBjb2RlIGNvbWluZyBmcm9tIHRoZSBjb2RlIHN0cmVhbS5cbiAgdmFyIGNvZGVfbWFzayA9ICgxIDw8IGN1cl9jb2RlX3NpemUpIC0gMTtcbiAgdmFyIGN1cl9zaGlmdCA9IDA7XG4gIHZhciBjdXIgPSAwO1xuXG4gIHZhciBvcCA9IDA7ICAvLyBPdXRwdXQgcG9pbnRlci5cbiAgXG4gIHZhciBzdWJibG9ja19zaXplID0gY29kZV9zdHJlYW1bcCsrXTtcblxuICAvLyBUT0RPKGRlYW5tKTogV291bGQgdXNpbmcgYSBUeXBlZEFycmF5IGJlIGFueSBmYXN0ZXI/ICBBdCBsZWFzdCBpdCB3b3VsZFxuICAvLyBzb2x2ZSB0aGUgZmFzdCBtb2RlIC8gYmFja2luZyBzdG9yZSB1bmNlcnRhaW50eS5cbiAgLy8gdmFyIGNvZGVfdGFibGUgPSBBcnJheSg0MDk2KTtcbiAgdmFyIGNvZGVfdGFibGUgPSBuZXcgSW50MzJBcnJheSg0MDk2KTsgIC8vIENhbiBiZSBzaWduZWQsIHdlIG9ubHkgdXNlIDIwIGJpdHMuXG5cbiAgdmFyIHByZXZfY29kZSA9IG51bGw7ICAvLyBUcmFjayBjb2RlLTEuXG5cbiAgd2hpbGUgKHRydWUpIHtcbiAgICAvLyBSZWFkIHVwIHRvIHR3byBieXRlcywgbWFraW5nIHN1cmUgd2UgYWx3YXlzIDEyLWJpdHMgZm9yIG1heCBzaXplZCBjb2RlLlxuICAgIHdoaWxlIChjdXJfc2hpZnQgPCAxNikge1xuICAgICAgaWYgKHN1YmJsb2NrX3NpemUgPT09IDApIGJyZWFrOyAgLy8gTm8gbW9yZSBkYXRhIHRvIGJlIHJlYWQuXG5cbiAgICAgIGN1ciB8PSBjb2RlX3N0cmVhbVtwKytdIDw8IGN1cl9zaGlmdDtcbiAgICAgIGN1cl9zaGlmdCArPSA4O1xuXG4gICAgICBpZiAoc3ViYmxvY2tfc2l6ZSA9PT0gMSkgeyAgLy8gTmV2ZXIgbGV0IGl0IGdldCB0byAwIHRvIGhvbGQgbG9naWMgYWJvdmUuXG4gICAgICAgIHN1YmJsb2NrX3NpemUgPSBjb2RlX3N0cmVhbVtwKytdOyAgLy8gTmV4dCBzdWJibG9jay5cbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIC0tc3ViYmxvY2tfc2l6ZTtcbiAgICAgIH1cbiAgICB9XG5cbiAgICAvLyBUT0RPKGRlYW5tKTogV2Ugc2hvdWxkIG5ldmVyIHJlYWxseSBnZXQgaGVyZSwgd2Ugc2hvdWxkIGhhdmUgcmVjZWl2ZWRcbiAgICAvLyBhbmQgRU9JLlxuICAgIGlmIChjdXJfc2hpZnQgPCBjdXJfY29kZV9zaXplKVxuICAgICAgYnJlYWs7XG5cbiAgICB2YXIgY29kZSA9IGN1ciAmIGNvZGVfbWFzaztcbiAgICBjdXIgPj49IGN1cl9jb2RlX3NpemU7XG4gICAgY3VyX3NoaWZ0IC09IGN1cl9jb2RlX3NpemU7XG5cbiAgICAvLyBUT0RPKGRlYW5tKTogTWF5YmUgc2hvdWxkIGNoZWNrIHRoYXQgdGhlIGZpcnN0IGNvZGUgd2FzIGEgY2xlYXIgY29kZSxcbiAgICAvLyBhdCBsZWFzdCB0aGlzIGlzIHdoYXQgeW91J3JlIHN1cHBvc2VkIHRvIGRvLiAgQnV0IGFjdHVhbGx5IG91ciBlbmNvZGVyXG4gICAgLy8gbm93IGRvZXNuJ3QgZW1pdCBhIGNsZWFyIGNvZGUgZmlyc3QgYW55d2F5LlxuICAgIGlmIChjb2RlID09PSBjbGVhcl9jb2RlKSB7XG4gICAgICAvLyBXZSBkb24ndCBhY3R1YWxseSBoYXZlIHRvIGNsZWFyIHRoZSB0YWJsZS4gIFRoaXMgY291bGQgYmUgYSBnb29kIGlkZWFcbiAgICAgIC8vIGZvciBncmVhdGVyIGVycm9yIGNoZWNraW5nLCBidXQgd2UgZG9uJ3QgcmVhbGx5IGRvIGFueSBhbnl3YXkuICBXZVxuICAgICAgLy8gd2lsbCBqdXN0IHRyYWNrIGl0IHdpdGggbmV4dF9jb2RlIGFuZCBvdmVyd3JpdGUgb2xkIGVudHJpZXMuXG5cbiAgICAgIG5leHRfY29kZSA9IGVvaV9jb2RlICsgMTtcbiAgICAgIGN1cl9jb2RlX3NpemUgPSBtaW5fY29kZV9zaXplICsgMTtcbiAgICAgIGNvZGVfbWFzayA9ICgxIDw8IGN1cl9jb2RlX3NpemUpIC0gMTtcblxuICAgICAgLy8gRG9uJ3QgdXBkYXRlIHByZXZfY29kZSA/XG4gICAgICBwcmV2X2NvZGUgPSBudWxsO1xuICAgICAgY29udGludWU7XG4gICAgfSBlbHNlIGlmIChjb2RlID09PSBlb2lfY29kZSkge1xuICAgICAgYnJlYWs7XG4gICAgfVxuXG4gICAgLy8gV2UgaGF2ZSBhIHNpbWlsYXIgc2l0dWF0aW9uIGFzIHRoZSBkZWNvZGVyLCB3aGVyZSB3ZSB3YW50IHRvIHN0b3JlXG4gICAgLy8gdmFyaWFibGUgbGVuZ3RoIGVudHJpZXMgKGNvZGUgdGFibGUgZW50cmllcyksIGJ1dCB3ZSB3YW50IHRvIGRvIGluIGFcbiAgICAvLyBmYXN0ZXIgbWFubmVyIHRoYW4gYW4gYXJyYXkgb2YgYXJyYXlzLiAgVGhlIGNvZGUgYmVsb3cgc3RvcmVzIHNvcnQgb2YgYVxuICAgIC8vIGxpbmtlZCBsaXN0IHdpdGhpbiB0aGUgY29kZSB0YWJsZSwgYW5kIHRoZW4gXCJjaGFzZXNcIiB0aHJvdWdoIGl0IHRvXG4gICAgLy8gY29uc3RydWN0IHRoZSBkaWN0aW9uYXJ5IGVudHJpZXMuICBXaGVuIGEgbmV3IGVudHJ5IGlzIGNyZWF0ZWQsIGp1c3QgdGhlXG4gICAgLy8gbGFzdCBieXRlIGlzIHN0b3JlZCwgYW5kIHRoZSByZXN0IChwcmVmaXgpIG9mIHRoZSBlbnRyeSBpcyBvbmx5XG4gICAgLy8gcmVmZXJlbmNlZCBieSBpdHMgdGFibGUgZW50cnkuICBUaGVuIHRoZSBjb2RlIGNoYXNlcyB0aHJvdWdoIHRoZVxuICAgIC8vIHByZWZpeGVzIHVudGlsIGl0IHJlYWNoZXMgYSBzaW5nbGUgYnl0ZSBjb2RlLiAgV2UgaGF2ZSB0byBjaGFzZSB0d2ljZSxcbiAgICAvLyBmaXJzdCB0byBjb21wdXRlIHRoZSBsZW5ndGgsIGFuZCB0aGVuIHRvIGFjdHVhbGx5IGNvcHkgdGhlIGRhdGEgdG8gdGhlXG4gICAgLy8gb3V0cHV0IChiYWNrd2FyZHMsIHNpbmNlIHdlIGtub3cgdGhlIGxlbmd0aCkuICBUaGUgYWx0ZXJuYXRpdmUgd291bGQgYmVcbiAgICAvLyBzdG9yaW5nIHNvbWV0aGluZyBpbiBhbiBpbnRlcm1lZGlhdGUgc3RhY2ssIGJ1dCB0aGF0IGRvZXNuJ3QgbWFrZSBhbnlcbiAgICAvLyBtb3JlIHNlbnNlLiAgSSBpbXBsZW1lbnRlZCBhbiBhcHByb2FjaCB3aGVyZSBpdCBhbHNvIHN0b3JlZCB0aGUgbGVuZ3RoXG4gICAgLy8gaW4gdGhlIGNvZGUgdGFibGUsIGFsdGhvdWdoIGl0J3MgYSBiaXQgdHJpY2t5IGJlY2F1c2UgeW91IHJ1biBvdXQgb2ZcbiAgICAvLyBiaXRzICgxMiArIDEyICsgOCksIGJ1dCBJIGRpZG4ndCBtZWFzdXJlIG11Y2ggaW1wcm92ZW1lbnRzICh0aGUgdGFibGVcbiAgICAvLyBlbnRyaWVzIGFyZSBnZW5lcmFsbHkgbm90IHRoZSBsb25nKS4gIEV2ZW4gd2hlbiBJIGNyZWF0ZWQgYmVuY2htYXJrcyBmb3JcbiAgICAvLyB2ZXJ5IGxvbmcgdGFibGUgZW50cmllcyB0aGUgY29tcGxleGl0eSBkaWQgbm90IHNlZW0gd29ydGggaXQuXG4gICAgLy8gVGhlIGNvZGUgdGFibGUgc3RvcmVzIHRoZSBwcmVmaXggZW50cnkgaW4gMTIgYml0cyBhbmQgdGhlbiB0aGUgc3VmZml4XG4gICAgLy8gYnl0ZSBpbiA4IGJpdHMsIHNvIGVhY2ggZW50cnkgaXMgMjAgYml0cy5cblxuICAgIHZhciBjaGFzZV9jb2RlID0gY29kZSA8IG5leHRfY29kZSA/IGNvZGUgOiBwcmV2X2NvZGU7XG5cbiAgICAvLyBDaGFzZSB3aGF0IHdlIHdpbGwgb3V0cHV0LCBlaXRoZXIge0NPREV9IG9yIHtDT0RFLTF9LlxuICAgIHZhciBjaGFzZV9sZW5ndGggPSAwO1xuICAgIHZhciBjaGFzZSA9IGNoYXNlX2NvZGU7XG4gICAgd2hpbGUgKGNoYXNlID4gY2xlYXJfY29kZSkge1xuICAgICAgY2hhc2UgPSBjb2RlX3RhYmxlW2NoYXNlXSA+PiA4O1xuICAgICAgKytjaGFzZV9sZW5ndGg7XG4gICAgfVxuXG4gICAgdmFyIGsgPSBjaGFzZTtcbiAgICBcbiAgICB2YXIgb3BfZW5kID0gb3AgKyBjaGFzZV9sZW5ndGggKyAoY2hhc2VfY29kZSAhPT0gY29kZSA/IDEgOiAwKTtcbiAgICBpZiAob3BfZW5kID4gb3V0cHV0X2xlbmd0aCkge1xuICAgICAgY29uc29sZS5sb2coXCJXYXJuaW5nLCBnaWYgc3RyZWFtIGxvbmdlciB0aGFuIGV4cGVjdGVkLlwiKTtcbiAgICAgIHJldHVybjtcbiAgICB9XG5cbiAgICAvLyBBbHJlYWR5IGhhdmUgdGhlIGZpcnN0IGJ5dGUgZnJvbSB0aGUgY2hhc2UsIG1pZ2h0IGFzIHdlbGwgd3JpdGUgaXQgZmFzdC5cbiAgICBvdXRwdXRbb3ArK10gPSBrO1xuXG4gICAgb3AgKz0gY2hhc2VfbGVuZ3RoO1xuICAgIHZhciBiID0gb3A7ICAvLyBUcmFjayBwb2ludGVyLCB3cml0aW5nIGJhY2t3YXJkcy5cblxuICAgIGlmIChjaGFzZV9jb2RlICE9PSBjb2RlKSAgLy8gVGhlIGNhc2Ugb2YgZW1pdHRpbmcge0NPREUtMX0gKyBrLlxuICAgICAgb3V0cHV0W29wKytdID0gaztcblxuICAgIGNoYXNlID0gY2hhc2VfY29kZTtcbiAgICB3aGlsZSAoY2hhc2VfbGVuZ3RoLS0pIHtcbiAgICAgIGNoYXNlID0gY29kZV90YWJsZVtjaGFzZV07XG4gICAgICBvdXRwdXRbLS1iXSA9IGNoYXNlICYgMHhmZjsgIC8vIFdyaXRlIGJhY2t3YXJkcy5cbiAgICAgIGNoYXNlID4+PSA4OyAgLy8gUHVsbCBkb3duIHRvIHRoZSBwcmVmaXggY29kZS5cbiAgICB9XG5cbiAgICBpZiAocHJldl9jb2RlICE9PSBudWxsICYmIG5leHRfY29kZSA8IDQwOTYpIHtcbiAgICAgIGNvZGVfdGFibGVbbmV4dF9jb2RlKytdID0gcHJldl9jb2RlIDw8IDggfCBrO1xuICAgICAgLy8gVE9ETyhkZWFubSk6IEZpZ3VyZSBvdXQgdGhpcyBjbGVhcmluZyB2cyBjb2RlIGdyb3d0aCBsb2dpYyBiZXR0ZXIuICBJXG4gICAgICAvLyBoYXZlIGFuIGZlZWxpbmcgdGhhdCBpdCBzaG91bGQganVzdCBoYXBwZW4gc29tZXdoZXJlIGVsc2UsIGZvciBub3cgaXRcbiAgICAgIC8vIGlzIGF3a3dhcmQgYmV0d2VlbiB3aGVuIHdlIGdyb3cgcGFzdCB0aGUgbWF4IGFuZCB0aGVuIGhpdCBhIGNsZWFyIGNvZGUuXG4gICAgICAvLyBGb3Igbm93IGp1c3QgY2hlY2sgaWYgd2UgaGl0IHRoZSBtYXggMTItYml0cyAodGhlbiBhIGNsZWFyIGNvZGUgc2hvdWxkXG4gICAgICAvLyBmb2xsb3csIGFsc28gb2YgY291cnNlIGVuY29kZWQgaW4gMTItYml0cykuXG4gICAgICBpZiAobmV4dF9jb2RlID49IGNvZGVfbWFzaysxICYmIGN1cl9jb2RlX3NpemUgPCAxMikge1xuICAgICAgICArK2N1cl9jb2RlX3NpemU7XG4gICAgICAgIGNvZGVfbWFzayA9IGNvZGVfbWFzayA8PCAxIHwgMTtcbiAgICAgIH1cbiAgICB9XG5cbiAgICBwcmV2X2NvZGUgPSBjb2RlO1xuICB9XG5cbiAgaWYgKG9wICE9PSBvdXRwdXRfbGVuZ3RoKSB7XG4gICAgY29uc29sZS5sb2coXCJXYXJuaW5nLCBnaWYgc3RyZWFtIHNob3J0ZXIgdGhhbiBleHBlY3RlZC5cIik7XG4gIH1cblxuICByZXR1cm4gb3V0cHV0O1xufVxuXG50cnkgeyBleHBvcnRzLkdpZldyaXRlciA9IEdpZldyaXRlcjsgZXhwb3J0cy5HaWZSZWFkZXIgPSBHaWZSZWFkZXIgfSBjYXRjaChlKSB7IH0gIC8vIENvbW1vbkpTLlxuIiwie0dpZlJlYWRlcn0gPSByZXF1aXJlICdvbWdnaWYnXG5cblBoYXNlci5Mb2FkZXIucHJvdG90eXBlLnNwcml0ZWdpZiA9IChrZXksIHVybCkgLT5cbiAgY29uc29sZS5sb2cgJ0dJRjogJywga2V5LCB1cmxcbiAgQGJpbmFyeSBrZXksIHVybCwgKFxuICAgIChrZXksIGRhdGEpID0+XG4gICAgICBzdGFydCA9IERhdGUubm93KClcbiAgICAgIHJlYWRlciA9IG5ldyBHaWZSZWFkZXIgbmV3IFVpbnQ4QXJyYXkgZGF0YVxuICAgICAgc3F1YXJlID0gTWF0aC5jZWlsIE1hdGguc3FydCByZWFkZXIubnVtRnJhbWVzKClcbiAgICAgIGNvbnNvbGUubG9nIHNxdWFyZSAqIHJlYWRlci53aWR0aCwgc3F1YXJlICogcmVhZGVyLmhlaWdodFxuICAgICAgY2FudmFzID0gUGhhc2VyLkNhbnZhcy5jcmVhdGUgc3F1YXJlICogcmVhZGVyLndpZHRoLCBzcXVhcmUgKiByZWFkZXIuaGVpZ2h0XG4gICAgICBjb250ZXh0ID0gY2FudmFzLmdldENvbnRleHQgJzJkJ1xuICAgICAgZm9yIGluZGV4IGluIFswLi4ucmVhZGVyLm51bUZyYW1lcygpXVxuICAgICAgICB4ID0gKGluZGV4ICUgc3F1YXJlKSAqIHJlYWRlci53aWR0aFxuICAgICAgICB5ID0gKE1hdGguZmxvb3IgaW5kZXggLyBzcXVhcmUpICogcmVhZGVyLmhlaWdodFxuICAgICAgICBpbWFnZSA9IG5ldyBJbWFnZURhdGEgcmVhZGVyLndpZHRoLCByZWFkZXIuaGVpZ2h0XG4gICAgICAgIGltYWdlW2tleV0gPSB2YWx1ZSBmb3Iga2V5LCB2YWx1ZSBvZiByZWFkZXIuZnJhbWVJbmZvIGluZGV4XG4gICAgICAgIGltYWdlLmRlbGF5ID0gaW1hZ2UuZGVsYXkgKiAxMCBpZiBpbWFnZS5kZWxheVxuICAgICAgICByZWFkZXIuZGVjb2RlQW5kQmxpdEZyYW1lUkdCQSBpbmRleCwgaW1hZ2UuZGF0YVxuICAgICAgICBjb250ZXh0LnB1dEltYWdlRGF0YSBpbWFnZSwgeCwgeVxuICAgICAgUGhhc2VyLkNhbnZhcy5hZGRUb0RPTSBjYW52YXMsIGRvY3VtZW50LmJvZHlcbiAgKSwgQGdhbWVcblxuIyB4aHIgPSBuZXcgWE1MSHR0cFJlcXVlc3RcbiMgeGhyLm9wZW4gJ0dFVCcsICdodHRwOi8vaS5pbWd1ci5jb20vOTBDczl0ci5naWYnXG4jIHhoci5yZXNwb25zZVR5cGUgPSAnYXJyYXlidWZmZXInXG4jIHhoci5vbmxvYWQgPSAoZXZlbnQpIC0+XG4jICAgdGltZSA9IERhdGUubm93KClcbiMgICByZWFkZXI9IG5ldyBHaWZSZWFkZXIgbmV3IFVpbnQ4QXJyYXkgQHJlc3BvbnNlXG4jICAgaW1hZ2VzPVxuIyAgICAgZm9yIGluZGV4IGluIFswLi4ucmVhZGVyLm51bUZyYW1lcygpXVxuIyAgICAgICBpbWFnZSA9IEBjcmVhdGVJbWFnZURhdGEgcmVhZGVyLndpZHRoLCByZWFkZXIuaGVpZ2h0XG4jICAgICAgIGltYWdlW2tleV0gPSB2YWx1ZSBmb3Iga2V5LCB2YWx1ZSBvZiByZWFkZXIuZnJhbWVJbmZvIGluZGV4XG4jICAgICAgIGltYWdlLmRlbGF5ID0gaW1hZ2UuZGVsYXkgKiAxMCBpZiBpbWFnZS5kZWxheVxuIyAgICAgICByZWFkZXIuZGVjb2RlQW5kQmxpdEZyYW1lUkdCQSBpbmRleCwgaW1hZ2UuZGF0YVxuIyAgICAgICBpbWFnZVxuIyAgIGltYWdlcy5sb29wQ291bnQ9IHJlYWRlci5sb29wQ291bnQoKVxuIyAgIGltYWdlcy5sb29wQ291bnQ/PSAtMVxuI1xuIyAgIGNvbnNvbGUubG9nIERhdGUubm93KCkgLSB0aW1lLCBpbWFnZXNcbiMgeGhyLnNlbmQoKVxuXG4jXG4jIHhociA9IG5ldyBYTUxIdHRwUmVxdWVzdFxuIyB4aHIub3BlbiAnR0VUJywgJ2h0dHA6Ly9pLmltZ3VyLmNvbS9uWnNrZVdULmdpZicsIHRydWVcbiMgIyBIYWNrIHRvIHBhc3MgYnl0ZXMgdGhyb3VnaCB1bnByb2Nlc3NlZC5cbiMgeGhyLm92ZXJyaWRlTWltZVR5cGUgJ3RleHQvcGxhaW47IGNoYXJzZXQ9eC11c2VyLWRlZmluZWQnXG4jXG4jIHhoci5vbnJlYWR5c3RhdGVjaGFuZ2UgPSAoZSkgLT5cbiMgICBpZiBAcmVhZHlTdGF0ZSA9PSA0IGFuZCBAc3RhdHVzID09IDIwMFxuIyAgICAgYmluU3RyID0gQHJlc3BvbnNlVGV4dFxuIyAgICAgaSA9IDBcbiMgICAgIGxlbiA9IGJpblN0ci5sZW5ndGhcbiMgICAgIHdoaWxlIGkgPCBsZW5cbiMgICAgICAgYyA9IGJpblN0ci5jaGFyQ29kZUF0KGkpXG4jICAgICAgICNTdHJpbmcuZnJvbUNoYXJDb2RlKGMgJiAweGZmKTtcbiMgICAgICAgYnl0ZSA9IGMgJiAweGZmXG4jICAgICAgICMgYnl0ZSBhdCBvZmZzZXQgaVxuIyAgICAgICArK2lcbiMgICByZXR1cm5cbiNcbiMgeGhyLnNlbmQoKVxuXG4jIC0tLVxuIyBnZW5lcmF0ZWQgYnkganMyY29mZmVlIDIuMC40XG5cbiMgZmlsZSA9XG4jICAgdHlwZTogJ2ltYWdlJ1xuIyAgIGtleTogJ2x1bmEnXG4jICAgdXJsOiAnaHR0cDovL2kuaW1ndXIuY29tL0xSUW5CaWgucG5nJ1xuIyAgIGRhdGE6IG51bGxcbiMgICBlcnJvcjogZmFsc2VcbiMgICBsb2FkZWQ6IGZhbHNlXG4jIGZpbGUuZGF0YSA9XG4jIGZpbGUuZGF0YS5uYW1lID0gZmlsZS5rZXlcbiNcbiMgZmlsZS5kYXRhLm9ubG9hZCA9IC0+XG4jXG4jICAgY29uc29sZS5sb2cgZmlsZS5kYXRhXG4jXG4jICAgZmlsZS5sb2FkZWQgPSB0cnVlXG4jICAgZ2FtZS5jYWNoZS5hZGRJbWFnZSBmaWxlLmtleSwgZmlsZS51cmwsIGZpbGUuZGF0YVxuIyAgIGx1bmEgPSBnYW1lLmFkZC5zcHJpdGUoZ2FtZS53b3JsZC5jZW50ZXJYLCBnYW1lLndvcmxkLmNlbnRlclksICdsdW5hJylcbiMgICBsdW5hLmFuY2hvci5zZXRUbyAwLjUsIDAuNVxuIyAgIHJldHVyblxuI1xuIyBmaWxlLmRhdGEub25lcnJvciA9IC0+XG4jICAgZmlsZS5lcnJvciA9IHRydWVcbiMgICByZXR1cm5cbiNcbiMgZmlsZS5kYXRhLmNyb3NzT3JpZ2luID0gJyonXG4jIGZpbGUuZGF0YS5zcmMgPSBmaWxlLnVybFxuIiwibHVuZSA9IHJlcXVpcmUgJ2x1bmUnXG5jb25zb2xlLmxvZyBsdW5lLnBoYXNlKClcblxucmVxdWlyZSAnLi9naWYuY29mZmVlJ1xuXG5jbGFzcyBMdW5hIGV4dGVuZHMgUGhhc2VyLkdhbWVcbiAgY29uc3RydWN0b3I6IC0+XG4gICAgdW5sZXNzIEAgaW5zdGFuY2VvZiBMdW5hIHRoZW4gcmV0dXJuIG5ldyBMdW5hXG4gICAgQGdhbWUgPSBzdXBlciB3aW5kb3cuaW5uZXJXaWR0aCwgd2luZG93LmlubmVySGVpZ2h0LCBQaGFzZXIuQVVUTywgJycsIEBcbiAgcHJlbG9hZDogLT5cbiAgICBAZ2FtZS5sb2FkLmNyb3NzT3JpZ2luID0gJyonXG4gICAgQGdhbWUubG9hZC5zcHJpdGVnaWYgJ204aFlFamInLCAnaHR0cDovL2kuaW1ndXIuY29tL204aFlFamIuZ2lmJ1xuXG4gIGNyZWF0ZTogLT5cblxud2luZG93LmdhbWUgPSBuZXcgTHVuYVxuIl19
