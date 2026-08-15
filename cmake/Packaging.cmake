include_guard(GLOBAL)

if(NOT TARGET bam-seek)
    message(FATAL_ERROR "BAM Seek packaging requires the bam-seek GUI target")
endif()

# App bundles belong at the root of a DMG. Windows executables and their DLLs
# also live together at the installer root so the loader can find them.
if(WIN32)
    install(TARGETS bam-seek
        RUNTIME_DEPENDENCIES
            DIRECTORIES
                "$<TARGET_FILE_DIR:bam-seek>"
                "$<TARGET_FILE_DIR:Qt6::Core>"
            PRE_EXCLUDE_REGEXES
                "api-ms-.*"
                "ext-ms-.*"
            POST_EXCLUDE_REGEXES
                ".*[\\\\/]System32[\\\\/].*\\.dll"
                ".*[\\\\/]SysWOW64[\\\\/].*\\.dll"
        RUNTIME DESTINATION . COMPONENT Runtime
        LIBRARY DESTINATION . COMPONENT Runtime
    )
else()
    install(TARGETS bam-seek
        BUNDLE DESTINATION . COMPONENT Runtime
        RUNTIME DESTINATION . COMPONENT Runtime
    )
endif()

# macdeployqt normally copies every image/input plugin in a Homebrew Qt
# installation. Some optional plugins can depend on Qt modules that this app
# neither uses nor requires. Seed only the native platform/style/network/TLS
# plugins, then ask macdeployqt to fix them up along with the app's libraries.
if(APPLE)
    qt_generate_deploy_script(
        TARGET bam-seek
        OUTPUT_SCRIPT bam_seek_deploy_script
        CONTENT [=[
set(bam_seek_bundle "$<TARGET_FILE_NAME:bam-seek>.app")
set(bam_seek_plugin_root
    "\${__QT_DEPLOY_QT_INSTALL_PREFIX}/\${__QT_DEPLOY_QT_INSTALL_PLUGINS}")
set(bam_seek_plugins
    "platforms/libqcocoa.dylib"
    "styles/libqmacstyle.dylib"
    "networkinformation/libqapplenetworkinformation.dylib"
    "tls/libqcertonlybackend.dylib"
    "tls/libqsecuretransportbackend.dylib"
)
foreach(bam_seek_plugin IN LISTS bam_seek_plugins)
    get_filename_component(
        bam_seek_plugin_source
        "\${bam_seek_plugin_root}/\${bam_seek_plugin}"
        REALPATH
    )
    if(NOT EXISTS "\${bam_seek_plugin_source}")
        message(FATAL_ERROR "Required Qt plugin not found: \${bam_seek_plugin_source}")
    endif()
    get_filename_component(bam_seek_plugin_type "\${bam_seek_plugin}" DIRECTORY)
    file(INSTALL
        DESTINATION
            "\${QT_DEPLOY_PREFIX}/\${bam_seek_bundle}/Contents/PlugIns/\${bam_seek_plugin_type}"
        TYPE FILE
        FILES "\${bam_seek_plugin_source}"
    )
endforeach()

qt_deploy_runtime_dependencies(
    EXECUTABLE "\${bam_seek_bundle}"
    NO_PLUGINS
    NO_TRANSLATIONS
    DEPLOY_TOOL_OPTIONS "-codesign=-"
)
]=]
    )
else()
    # windeployqt copies Qt DLLs plus the platform, image, network, and TLS
    # plugins selected for the application.
    qt_generate_deploy_app_script(
        TARGET bam-seek
        OUTPUT_SCRIPT bam_seek_deploy_script
        NO_TRANSLATIONS
        NO_UNSUPPORTED_PLATFORM_ERROR
    )
endif()
install(SCRIPT "${bam_seek_deploy_script}" COMPONENT Runtime)

set(CPACK_PACKAGE_NAME "BAM Seek")
set(CPACK_PACKAGE_VENDOR "BAM Seek")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "Desktop application for reviewing variant allele frequencies in BAM files")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/sa501428/bam-seek")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "BAM Seek")
set(CPACK_PACKAGE_DIRECTORY "${CMAKE_BINARY_DIR}/packages")
set(CPACK_PACKAGE_CHECKSUM "SHA256")
set(CPACK_MONOLITHIC_INSTALL ON)
set(CPACK_STRIP_FILES OFF)

if(CMAKE_OSX_ARCHITECTURES)
    string(REPLACE ";" "-" bam_seek_package_arch "${CMAKE_OSX_ARCHITECTURES}")
elseif(CMAKE_SYSTEM_PROCESSOR)
    set(bam_seek_package_arch "${CMAKE_SYSTEM_PROCESSOR}")
else()
    set(bam_seek_package_arch "unknown")
endif()

if(APPLE)
    set(CPACK_GENERATOR "DragNDrop")
    set(CPACK_PACKAGE_FILE_NAME
        "BAM-Seek-${PROJECT_VERSION}-macOS-${bam_seek_package_arch}")
    set(CPACK_DMG_VOLUME_NAME "BAM Seek ${PROJECT_VERSION}")
    set(CPACK_DMG_FORMAT "UDZO")
elseif(WIN32)
    set(BAM_SEEK_WINDOWS_GENERATORS "NSIS;WIX" CACHE STRING
        "Windows CPack generators (NSIS for .exe and WIX for .msi)")
    set(CPACK_GENERATOR "${BAM_SEEK_WINDOWS_GENERATORS}")
    set(CPACK_PACKAGE_FILE_NAME
        "BAM-Seek-${PROJECT_VERSION}-Windows-x64")

    # Both generators use this list to create Start Menu launchers. NSIS also
    # creates an optional desktop shortcut through CPACK_CREATE_DESKTOP_LINKS.
    set(CPACK_PACKAGE_EXECUTABLES "bam-seek" "BAM Seek")
    set(CPACK_CREATE_DESKTOP_LINKS "bam-seek")

    set(CPACK_NSIS_DISPLAY_NAME "BAM Seek ${PROJECT_VERSION}")
    set(CPACK_NSIS_PACKAGE_NAME "BAM Seek")
    set(CPACK_NSIS_INSTALLED_ICON_NAME "bam-seek.exe")
    set(CPACK_NSIS_MUI_FINISHPAGE_RUN "bam-seek.exe")
    set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
    set(CPACK_NSIS_MANIFEST_DPI_AWARE ON)
    set(CPACK_NSIS_UNINSTALL_NAME "Uninstall BAM Seek")
    set(CPACK_NSIS_HELP_LINK "${CPACK_PACKAGE_HOMEPAGE_URL}")
    set(CPACK_NSIS_URL_INFO_ABOUT "${CPACK_PACKAGE_HOMEPAGE_URL}")

    # Keep this GUID stable across releases so Windows Installer upgrades an
    # existing BAM Seek installation instead of installing a second product.
    set(CPACK_WIX_UPGRADE_GUID "C09CB1C5-B734-4090-8A10-0DDD6470AAC7")
    set(CPACK_WIX_PROGRAM_MENU_FOLDER "BAM Seek")
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
        set(CPACK_WIX_VERSION 4)
    endif()
else()
    message(WARNING "Native BAM Seek packages are only configured for macOS and Windows")
endif()

include(CPack)
