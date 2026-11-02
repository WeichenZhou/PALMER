TARGET := PALMER
SRC := $(TARGET).cpp
PIPELINE_SRC := $(wildcard scp/*.cpp scp/*.hpp)
PKG_CONFIG ?= pkg-config
HTSLIB_MIN_VERSION := 1.24
HTS_CFLAGS = $(shell $(PKG_CONFIG) --cflags htslib 2>/dev/null)
HTS_LIBS = $(shell $(PKG_CONFIG) --libs htslib 2>/dev/null)

$(TARGET): $(SRC) $(PIPELINE_SRC) makefile | check-htslib
	$(CXX) -o $@ $(HTS_CFLAGS) $(CXXFLAGS) $< -O3 -w -std=c++17 -pthread -lstdc++fs $(HTS_LIBS)

check-htslib:
	@if ! command -v $(PKG_CONFIG) >/dev/null 2>&1; then \
		echo "ERROR: pkg-config is required to locate HTSlib $(HTSLIB_MIN_VERSION) or newer."; \
		exit 1; \
	fi
	@if ! $(PKG_CONFIG) --exists htslib; then \
		echo "ERROR: HTSlib $(HTSLIB_MIN_VERSION) or newer was not found. Install its development files and set PKG_CONFIG_PATH if needed."; \
		exit 1; \
	fi
	@if ! $(PKG_CONFIG) --atleast-version=$(HTSLIB_MIN_VERSION) htslib; then \
		echo "ERROR: HTSlib $(HTSLIB_MIN_VERSION) or newer is required; found $$($(PKG_CONFIG) --modversion htslib)."; \
		exit 1; \
	fi
	@echo "Using HTSlib $$($(PKG_CONFIG) --modversion htslib)"

.PHONY: check-htslib clean
clean:
	rm -f $(TARGET)
