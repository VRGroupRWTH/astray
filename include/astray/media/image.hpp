#pragma once

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <vector>

#include <astray/math/indexing.hpp>
#include <astray/math/linear_algebra.hpp>
#include <astray/utility/warning.hpp>

DISABLE_WARNING_PUSH
DISABLE_WARNING_PRAGMAS
DISABLE_WARNING_NO_DEPRECATED
#include <astray/third_party/stb/stb_image.h>
#include <astray/third_party/stb/stb_image_write.h>
DISABLE_WARNING_POP

namespace ast
{
template <typename type>
class image
{
public:
  using storage_type = std::vector<type>;
  using size_type    = vector2<std::int32_t>;
  
  explicit image  (const size_type& size = {32, 32}, const type& value = type::Constant(255)) : data(size.prod(), value), size(size)
  {
    
  }
  explicit image  (const std::filesystem::path& filepath)
  {
    load(filepath);
  }
  image           (const image&       that) = default;
  image           (      image&&      temp) = default;
 ~image           ()                        = default;
  image& operator=(const image&       that) = default;
  image& operator=(      image&&      temp) = default;
  
  constexpr       type& at(const size_type& index)
  {
    return data[ravel_multi_index<size_type, true>(index, size)];
  }
  constexpr const type& at(const size_type& index) const
  {
    return data[ravel_multi_index<size_type, true>(index, size)];
  }

  void load(const std::filesystem::path& filepath)
  {
    if (!exists(filepath))
      throw std::runtime_error("File does not exist!");

    auto components = 0;
    auto* const raw = reinterpret_cast<type*>(stbi_load(filepath.string().c_str(), &size[0], &size[1], &components, sizeof(type)));
    data = std::vector<type>(raw, raw + size[0] * size[1]); // Copy will not be necessary if https://github.com/nothings/stb/issues/58 is resolved.
    stbi_image_free(raw);
  }
  void save(const std::filesystem::path& filepath) const
  {
    const auto extension = filepath.extension();
    if      (extension == ".bmp") stbi_write_bmp(filepath.string().c_str(), size[0], size[1], sizeof(type), data.data());
    else if (extension == ".jpg") stbi_write_jpg(filepath.string().c_str(), size[0], size[1], sizeof(type), data.data(), 100);
    else if (extension == ".png") stbi_write_png(filepath.string().c_str(), size[0], size[1], sizeof(type), data.data(), size[0] * sizeof(type));
    else if (extension == ".tga") stbi_write_tga(filepath.string().c_str(), size[0], size[1], sizeof(type), data.data());
    else throw std::runtime_error("Unsupported extension!");
  }
  
  std::vector<type>     data;
  vector2<std::int32_t> size;
};
}
