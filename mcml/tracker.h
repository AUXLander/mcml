#pragma once
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <shared_mutex>
#include <ios>
#include <cassert>
#include <map>
#include <set>
#include <unordered_map>

#include <boost/functional/hash.hpp>


template<class T>
struct point
{
	T x;
	T y;
	T z;

	size_t layer;

	point(T x, T y, T z, size_t layer) :
		x(x), y(y), z(z), layer(layer)
	{;}

	point(const point<T>& other) :
		x(other.x), y(other.y), z(other.z),
		layer(other.layer)
	{;}

	friend std::fstream& operator<<(std::fstream& stream, const point<T>& p)
	{
		stream.write((const char*)&p.x, sizeof(x));
		stream.write((const char*)&p.y, sizeof(y));
		stream.write((const char*)&p.z, sizeof(z));

		stream.write((const char*)&p.layer, sizeof(layer));

		return stream;
	}

	bool operator<(const point<T>& other) const
	{
		return layer < other.layer;
	}
};

using point_t = point<double>;
//using unique_point_t = std::unique_ptr<>;


constexpr static size_t dpi_x = 100U;
constexpr static size_t dpi_y = 100U;
constexpr static size_t dpi_z = 100U;


constexpr static size_t dots_per_x = 1U;
constexpr static size_t dots_per_y = dpi_x * dots_per_x;
constexpr static size_t dots_per_z = dpi_y * dots_per_y;
constexpr static size_t dots_per_l = dpi_z * dots_per_z;

struct p_point
{
	size_t x_idx;
	size_t y_idx;
	size_t z_idx;
	size_t l_idx;

	p_point(size_t x_idx, size_t y_idx, size_t z_idx, size_t l_idx) :
		x_idx(x_idx), y_idx(y_idx), z_idx(z_idx), l_idx(l_idx)
	{;}

	size_t to_index() const
	{
		return dots_per_x * x_idx +
			   dots_per_y * y_idx +
			   dots_per_z * z_idx + 
			   dots_per_l * l_idx;
	}

	bool operator==(const p_point& other) const
	{
		return x_idx == other.x_idx && 
			   y_idx == other.y_idx && 
			   z_idx == other.z_idx && 
			   l_idx == other.l_idx;
	}

	struct hash
	{
		size_t operator()(const p_point& p) const
		{
			size_t hash = 0;

			boost::hash_combine(hash, p.x_idx);
			boost::hash_combine(hash, p.y_idx);
			boost::hash_combine(hash, p.z_idx);
			boost::hash_combine(hash, p.l_idx);

			return hash;
		}
	};
};




struct tracker
{
	class id_t
	{
		size_t __val{0U};

	public:
		size_t next()
		{
			return __val++;
		}
	};

	struct local_thread_storage
	{
		using dot_counter_t = size_t;

		constexpr static double min_x = -1.0;
		constexpr static double min_y = -1.0;
		constexpr static double min_z = -1.0;

		constexpr static double max_x = +1.0;
		constexpr static double max_y = +1.0;
		constexpr static double max_z = +2.0;

		constexpr static double step_x = (max_x - min_x) / static_cast<double>(dpi_x);
		constexpr static double step_y = (max_y - min_y) / static_cast<double>(dpi_y);
		constexpr static double step_z = (max_z - min_z) / static_cast<double>(dpi_z);

		using dot_storage_t = std::vector<dot_counter_t>;

		friend class tracker;

		enum class combine_mode_e
		{
			COMBINE_ADD,
		};

	private:

		dot_storage_t __heat_map;

		std::unordered_map<p_point, size_t*, p_point::hash> __heat_hash;

		size_t __total_count;

	public:
		local_thread_storage() :
			__heat_map{},
			__total_count {0}
		{
			// 6 is count of layers
			__heat_map.resize(dpi_x* dpi_y* dpi_z* 6U, 0U);
		}

		local_thread_storage(const local_thread_storage& other) = delete;

		local_thread_storage(local_thread_storage&& other) = default;

		void track(const double& x, const double& y, const double& z, size_t layer)
		{
			const p_point key {
				static_cast<dot_counter_t>((x - min_x) / step_x),
				static_cast<dot_counter_t>((y - min_y) / step_y),
				static_cast<dot_counter_t>((z - min_z) / step_z),
				layer
			};

			if (key.to_index() > __heat_map.size())
			{
				key.to_index();
			}

			auto& item = __heat_map[key.to_index()];

			__heat_hash.try_emplace(key, &item);

			++item;
			++__total_count;
		}

		void combine(const local_thread_storage& lhs, combine_mode_e mode)
		{
			assert(__heat_map.size() == lhs.__heat_map.size());

			switch (mode)
			{
				case combine_mode_e::COMBINE_ADD:
				{
					auto coefficient = (float)lhs.__heat_hash.size() / (float)lhs.__heat_map.size();

					if (coefficient <= 0.25)
					{
						for (auto& [key, value] : lhs.__heat_hash)
						{
							auto& item = __heat_map[key.to_index()];

							__heat_hash.try_emplace(key, &item);

							item += *value;
						}
					}
					else
					{
						intptr_t size = __heat_map.size();

						#pragma omp parallel for
						for (intptr_t index = 0; index < size; ++index)
						{
							__heat_map[index] += lhs.__heat_map[index];
						}
					}

					__total_count += lhs.__total_count;
				};
				break;
			}
		}

		void write(std::fstream& stream) const 
		{
			stream.write((const char*)&__total_count, sizeof(__total_count));

			for (const auto& point : __heat_map)
			{
				stream << point;
			}
		}
	};

private:

	local_thread_storage __main_thread_storage;

	constexpr static auto openmode = std::ios_base::binary | std::ios_base::out;

	struct fstream_deleter
	{
		void operator()(std::fstream* s)
		{
			assert(s);
			s->close();
			delete s;
		}
	};

	std::unique_ptr<std::fstream, fstream_deleter> __file;

	mutable std::shared_mutex __mutex;

	tracker() {;}
	
public:

	static tracker& instance()
	{
		static tracker INSTANCE;
		return INSTANCE;
	}

	void track(local_thread_storage&& data)
	{
		std::unique_lock lock(__mutex);

		__main_thread_storage.combine(std::move(data), local_thread_storage::combine_mode_e::COMBINE_ADD);
	}

	void set_file(const char* filename)
	{
		__file.reset(new std::fstream{ filename, openmode });
	}

	template<typename Tfunc>
	void set_headers(Tfunc fgen)
	{
		std::shared_lock lock(__mutex);

		assert(__file);

		if (__file)
		{
			fgen(*__file);

			__file->flush();
		}
	}

	void write() const
	{
		std::shared_lock lock(__mutex);

		assert(__file);

		if (__file)
		{
			id_t idgen;
			size_t nums = 0U;

			double min_x = -1.0;
			double min_y = -1.0;
			double min_z = -1.0;

			double max_x = +1.0;
			double max_y = +1.0;
			double max_z = +2.0;

			auto& stream = *__file;

			__file->write((const char*)&__main_thread_storage.__total_count, sizeof(__main_thread_storage.__total_count));

			__file->flush();

			__main_thread_storage.write(stream);

			__file->flush();

			__file->seekg(sizeof(short) + sizeof(long));
			__file->seekp(sizeof(short) + sizeof(long));

			__file->write((const char*)&dpi_x, sizeof(dpi_x));
			__file->write((const char*)&dpi_y, sizeof(dpi_y));
			__file->write((const char*)&dpi_z, sizeof(dpi_z));

			__file->write((const char*)&min_x, sizeof(min_x));
			__file->write((const char*)&min_y, sizeof(min_y));
			__file->write((const char*)&min_z, sizeof(min_z));

			__file->write((const char*)&max_x, sizeof(max_x));
			__file->write((const char*)&max_y, sizeof(max_y));
			__file->write((const char*)&max_z, sizeof(max_z));

			__file->flush();
		}
	}
};

