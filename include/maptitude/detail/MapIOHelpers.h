#ifndef MAPTITUDE_DETAIL_MAPIOHELPERS_H
#define MAPTITUDE_DETAIL_MAPIOHELPERS_H

#include <filesystem>
#include <functional>

namespace Maptitude {
namespace detail {

/// Draws reserve_temporary_sibling makes before giving up.
constexpr int TEMPORARY_NAME_ATTEMPTS = 8;

/// Reserve a fresh file beside @p dest by exclusive create and return its path.
///
/// The name is ".maptitude-<pid>-<hex><extension of dest>", with the hex drawn
/// from @p entropy once per attempt. A candidate that already exists -- a plain
/// file, or a symbolic link whatever it resolves to -- is skipped for a fresh
/// draw rather than truncated, up to TEMPORARY_NAME_ATTEMPTS draws. write_map
/// draws from std::random_device; the parameter exists so that a collision on
/// an exact candidate can be produced under test, which 32 bits of device
/// entropy otherwise rule out.
///
/// @throws GridError naming @p dest if the create fails for a reason other
///         than an existing entry, or once every draw has collided.
std::filesystem::path reserve_temporary_sibling(
    const std::filesystem::path& dest,
    const std::function<unsigned int()>& entropy);

}  // namespace detail
}  // namespace Maptitude

#endif  // MAPTITUDE_DETAIL_MAPIOHELPERS_H
