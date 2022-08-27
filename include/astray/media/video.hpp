#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>

#ifdef ASTRAY_USE_FFMPEG
extern "C"
{
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/opt.h>
#include <libswscale/swscale.h>
}
#endif

#include <astray/math/linear_algebra.hpp>
#include <astray/media/image.hpp>
#include <astray/utility/warning.hpp>

DISABLE_WARNING_PUSH
DISABLE_WARNING_PRAGMAS
DISABLE_WARNING_NO_DEPRECATED

namespace ast
{
class video
{
public:
  explicit video  (
    const std::filesystem::path&      filepath   , 
    const vector2<std::int32_t>&      size       = vector2<std::int32_t>(1920, 1080),
    const std::int32_t                frame_rate = 30, 
    const std::int64_t                bit_rate   = 10000000,
    const std::optional<std::string>& codec_name = std::nullopt)
  {
#ifdef ASTRAY_USE_FFMPEG
    av_register_all();

    output_format_ = av_guess_format(nullptr, filepath.string().c_str(), nullptr);
    if (!output_format_) 
      throw std::runtime_error("Failed av_guess_format.");

    if(avformat_alloc_output_context2(&format_context_, output_format_, nullptr, nullptr) < 0)
      throw std::runtime_error("Failed avformat_alloc_output_context2.");

    auto* codec = avcodec_find_encoder(output_format_->video_codec);
    if (codec_name)
      codec = avcodec_find_encoder_by_name(codec_name.value().c_str());
    if (!codec)
      throw std::runtime_error("Failed avcodec_find_encoder.");

    stream_ = avformat_new_stream(format_context_, codec);
    if (!stream_)
      throw std::runtime_error("Failed avformat_new_stream.");

    codec_context_            = stream_->codec;
    codec_context_->width     = size[0];
    codec_context_->height    = size[1];
    codec_context_->pix_fmt   = AV_PIX_FMT_YUV420P;
    codec_context_->time_base = AVRational {1, frame_rate};
    codec_context_->bit_rate  = bit_rate;
    if (format_context_->oformat->flags & AVFMT_GLOBALHEADER)
      codec_context_->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

    AVDictionary* codec_options = nullptr;
    av_dict_set(&codec_options, "preset", "slow", 0);
    av_dict_set(&codec_options, "crf"   , "20"  , 0);
    if (avcodec_open2(codec_context_, codec, &codec_options) < 0)
      throw std::runtime_error("Failed avcodec_open2.");
    av_dict_free(&codec_options);

    stream_->time_base = codec_context_->time_base;
    av_dump_format(format_context_, 0, filepath.string().c_str(), 1);
    if (!(output_format_->flags & AVFMT_NOFILE) && avio_open(&format_context_->pb, filepath.string().c_str(), AVIO_FLAG_WRITE) < 0)
      throw std::runtime_error("Failed avio_open.");

    if (avformat_write_header(format_context_, nullptr) < 0)
      throw std::runtime_error("Failed avformat_write_header.");

    frame_         = av_frame_alloc();
    frame_->format = AV_PIX_FMT_YUV420P;
    frame_->width  = codec_context_->width;
    frame_->height = codec_context_->height;
    if (av_frame_get_buffer(frame_, 0) < 0)
      throw std::runtime_error("Failed av_frame_get_buffer.");

    sws_context_ = sws_getContext(size[0], size[1], AV_PIX_FMT_RGB24, size[0], size[1], AV_PIX_FMT_YUV420P, SWS_FAST_BILINEAR, nullptr, nullptr, nullptr);
    if (!sws_context_)
      throw std::runtime_error("Failed sws_getContext.");
#endif
  }
  video           (const video&  that) = delete ;
  video           (      video&& temp) = default;
  virtual ~video  ()
  {
#ifdef ASTRAY_USE_FFMPEG
    // Receive delayed frames.
    AVPacket packet;
    av_init_packet(&packet);
    packet.data = nullptr;
    packet.size = 0;

    auto received = 1;
    while (received)
    {
      avcodec_encode_video2(codec_context_, &packet, nullptr, &received);
      if (received)
      {
        std::fflush(stdout);
        av_packet_rescale_ts(&packet, codec_context_->time_base, stream_->time_base);
        packet.stream_index = stream_->index;
        av_interleaved_write_frame(format_context_, &packet);
        av_packet_unref(&packet);
      }
    }

    av_write_trailer(format_context_);

    if (!(output_format_->flags & AVFMT_NOFILE))
      avio_closep(&format_context_->pb);
    avcodec_close(stream_->codec);

    sws_freeContext      (sws_context_   );
    avformat_free_context(format_context_);
    av_frame_free        (&frame_        );
#endif
  }
  video& operator=(const video&  that) = delete ;
  video& operator=(      video&& temp) = default;

  void append(image<vector3<std::uint8_t>>& frame)
  {
#ifdef ASTRAY_USE_FFMPEG
    std::uint8_t*      data     [8] {&frame.data[0][0]};
    const std::int32_t line_size[8] {static_cast<std::int32_t>(frame.size[0] * sizeof(vector3<std::uint8_t>))};
    sws_scale(sws_context_, data, line_size, 0, frame.size[1], frame_->data, frame_->linesize);

    AVPacket packet;
    av_init_packet(&packet);
    packet.data = nullptr;
    packet.size = 0;

    frame_->pts = frame_count_++;

    auto received = 1;
    if (avcodec_encode_video2(codec_context_, &packet, frame_, &received) < 0)
      throw std::runtime_error("Failed avcodec_encode_video2.");
    if (received)
    {
      std::fflush(stdout);
      av_packet_rescale_ts(&packet, codec_context_->time_base, stream_->time_base);
      packet.stream_index = stream_->index;

      if (av_interleaved_write_frame(format_context_, &packet) < 0)
        throw std::runtime_error("Failed av_interleaved_write_frame.");
      av_packet_unref(&packet);
    }
#endif
  }

protected:
#ifdef ASTRAY_USE_FFMPEG
  AVOutputFormat*  output_format_  = nullptr;
  AVFormatContext* format_context_ = nullptr;
  AVCodecContext*  codec_context_  = nullptr;
  AVStream*        stream_         = nullptr;
  AVFrame*         frame_          = nullptr;
  SwsContext*      sws_context_    = nullptr;
#endif
  std::int64_t     frame_count_    = 0;
};
}

DISABLE_WARNING_POP
