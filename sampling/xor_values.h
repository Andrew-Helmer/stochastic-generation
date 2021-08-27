/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * These are xor_values for the various sequences we implement.
 */
#ifndef SAMPLING_XOR_VALUES_H
#define SAMPLING_XOR_VALUES_H

#include <cstdint>

namespace sampling {

static constexpr uint32_t pmj02_xors[2][32] = {
    {0x0, 0x0, 0x2, 0x6, 0x6, 0xe, 0x36, 0x4e, 0x16, 0x2e, 0x276, 0x6ce, 0x716, 0xc2e, 0x3076, 0x40ce, 0x116, 0x22e, 0x20676, 0x60ece, 0x61716, 0xe2c2e, 0x367076, 0x4ec0ce, 0x170116, 0x2c022e, 0x2700676, 0x6c00ece, 0x7001716, 0xc002c2e, 0x30007076, 0x4000c0ce},
    {0x0, 0x1, 0x3, 0x3, 0x7, 0x1b, 0x27, 0xb, 0x17, 0x13b, 0x367, 0x38b, 0x617, 0x183b, 0x2067, 0x8b, 0x117, 0x1033b, 0x30767, 0x30b8b, 0x71617, 0x1b383b, 0x276067, 0xb808b, 0x160117, 0x138033b, 0x3600767, 0x3800b8b, 0x6001617, 0x1800383b, 0x20006067, 0x808b}
};

/*
 * Derived from the direction numbers of Joe & Kuo (2008), generated from the
 * generator matrices provided with the PBRT renderer.
 */
static constexpr int MAX_SOBOL_DIM = 64;
static constexpr uint32_t sobol_xors[MAX_SOBOL_DIM][32] = {
    {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},
    {0x0, 0x1, 0x1, 0x7, 0x1, 0x13, 0x15, 0x7f, 0x1, 0x103, 0x105, 0x70f, 0x111, 0x1333, 0x1555, 0x7fff, 0x1, 0x10003, 0x10005, 0x7000f, 0x10011, 0x130033, 0x150055, 0x7f00ff, 0x10101, 0x1030303, 0x1050505, 0x70f0f0f, 0x1111111, 0x13333333, 0x15555555, 0x7fffffff},
    {0x0, 0x1, 0x3, 0x1, 0x5, 0x1f, 0x2b, 0x3d, 0x11, 0x133, 0x377, 0x199, 0x445, 0x1ccf, 0x2ddb, 0x366d, 0x101, 0x10303, 0x30707, 0x10909, 0x51515, 0x1f3f3f, 0x2b6b6b, 0x3dbdbd, 0x101011, 0x1303033, 0x3707077, 0x1909099, 0x4515145, 0x1cf3f3cf, 0x2db6b6db, 0x36dbdb6d},
    {0x0, 0x1, 0x0, 0x3, 0xd, 0xc, 0x5, 0x4f, 0x14, 0xe7, 0x329, 0x39c, 0x11, 0x1033, 0x44, 0x30bb, 0xd1cd, 0xc2ec, 0x5415, 0x4fc3f, 0x15054, 0xe5c97, 0x32e5b9, 0x39725c, 0x101, 0x1000303, 0x404, 0x3000b0b, 0xd001d1d, 0xc002c2c, 0x5004545, 0x4f00cfcf},
    {0x0, 0x0, 0x0, 0x5, 0xa, 0x14, 0x11, 0x22, 0x44, 0x19d, 0x33a, 0x674, 0x101, 0x202, 0x404, 0x5d0d, 0xba1a, 0x17434, 0x14151, 0x282a2, 0x50544, 0x1a4e9d, 0x349d3a, 0x693a74, 0x10001, 0x20002, 0x40004, 0x50d000d, 0xa1a001a, 0x14340034, 0x11510051, 0x22a200a2},
    {0x0, 0x0, 0x2, 0x6, 0x3, 0x6, 0x2a, 0x72, 0x5, 0xa, 0x21e, 0x636, 0x35f, 0x6be, 0x2bc2, 0x713a, 0x11, 0x22, 0x20066, 0x600ee, 0x30123, 0x60246, 0x2a06ca, 0x720fd2, 0x51155, 0xa22aa, 0x21e67fe, 0x636ed56, 0x35e26af, 0x6bc4d5e, 0x2bc4d7e2, 0x7135e29a},
    {0x0, 0x1, 0x1, 0x2, 0x9, 0xb, 0x3d, 0x7a, 0x41, 0x1c3, 0x45, 0x8a, 0xf59, 0x1eb, 0x223d, 0x447a, 0x1001, 0x13003, 0x15005, 0x2a00a, 0x89019, 0x9b02b, 0x3ad07d, 0x75a0fa, 0x551141, 0x1ff33c3, 0x15445, 0x2a88a, 0xeac8f59, 0x3f591eb, 0x241eb23d, 0x483d647a},
    {0x0, 0x0, 0x1, 0x2, 0x1, 0x5, 0xa, 0x31, 0x62, 0x75, 0x11, 0x22, 0x455, 0x8aa, 0x501, 0x1675, 0x2cea, 0xcfa1, 0x19f42, 0x1f125, 0x101, 0x202, 0x100505, 0x200a0a, 0x101111, 0x502525, 0xa04a4a, 0x310b1b1, 0x6216362, 0x7527775, 0x1141511, 0x2282a22},
    {0x0, 0x0, 0x1, 0x2, 0x5, 0x9, 0x12, 0xd, 0x1a, 0x1d, 0x41, 0x82, 0x545, 0xa8a, 0x1155, 0x2e69, 0x5cd2, 0x17cd, 0x2f9a, 0xf15d, 0x1001, 0x2002, 0x105005, 0x20a00a, 0x515015, 0x929029, 0x1252052, 0xd8d08d, 0x1b1a11a, 0x1f1d21d, 0x4541441, 0x8a82882},
    {0x0, 0x0, 0x3, 0x6, 0xf, 0xf, 0x1e, 0x4d, 0x9a, 0x145, 0x55, 0xaa, 0xdab, 0x1b56, 0x3a53, 0x35a3, 0x6b46, 0x10869, 0x210d2, 0x5ff41, 0x1111, 0x2222, 0x307777, 0x60eeee, 0xf1eeef, 0xf2dddf, 0x1e5bbbe, 0x4dc111d, 0x9b8223a, 0x14672215, 0x5114405, 0xa22880a},
    {0x0, 0x0, 0x1, 0x0, 0x0, 0x17, 0x2e, 0x6b, 0xb8, 0x170, 0x115, 0x22a, 0x141, 0x8a8, 0x1150, 0x689b, 0xd136, 0x14af7, 0x344d8, 0x689b0, 0x10111, 0x20222, 0x150555, 0x80888, 0x101110, 0x1473447, 0x28e688e, 0x65be55b, 0xa39a238, 0x14734470, 0x14404145, 0x2880828a},
    {0x0, 0x0, 0x0, 0x4, 0xe, 0x1b, 0x36, 0x6c, 0x34, 0xf2, 0x145, 0x28a, 0x514, 0x1f3c, 0x21e6, 0x5917, 0xb22e, 0x1645c, 0x1ace4, 0x18fba, 0x11011, 0x22022, 0x44044, 0x4cc0cc, 0xffe1fe, 0x188b38b, 0x3116716, 0x622ce2c, 0x2675274, 0xdfd0dd2, 0x11410115, 0x2282022a},
    {0x0, 0x1, 0x1, 0x3, 0xc, 0x1d, 0x7, 0x49, 0xaf, 0xcc, 0x151, 0x7f3, 0x15, 0x57b, 0x2adc, 0x45ad, 0x4ef7, 0x1d319, 0x3639f, 0x1610c, 0x11101, 0x133303, 0x155505, 0x3bbb0b, 0xdddc1c, 0x1eeed3d, 0x333747, 0x45559c9, 0xb445eaf, 0xe220ecc, 0x10114451, 0x7033ccf3},
    {0x0, 0x1, 0x3, 0x1, 0xe, 0x2, 0x3, 0x45, 0xc9, 0x5b, 0x3a2, 0xe6, 0x5, 0x100f, 0x301b, 0x102d, 0xe066, 0x20aa, 0x314f, 0x453d1, 0xc96ed, 0x5bb37, 0x3a392a, 0xe4b7e, 0x11, 0x1000033, 0x3000077, 0x1000099, 0xe0001fe, 0x2000202, 0x3000473, 0x45000c95},
    {0x0, 0x0, 0x0, 0x7, 0x5, 0xd, 0x1b, 0x36, 0x6c, 0x159, 0x87, 0x8f, 0x145, 0x28a, 0x514, 0x7cf3, 0x4001, 0xf6d9, 0x1fcf7, 0x3f9ee, 0x7f3dc, 0x16117d, 0xdc05b, 0x27673, 0x11011, 0x22022, 0x44044, 0x70ff0ff, 0x5145145, 0xd2fd2fd, 0x1b5eb5eb, 0x36bd6bd6},
    {0x0, 0x1, 0x0, 0x2, 0x8, 0x2, 0x21, 0x23, 0x84, 0x1ca, 0x118, 0x4e2, 0x401, 0x1c03, 0x1004, 0x80a, 0xe018, 0xa822, 0x39461, 0xbca3, 0xe5184, 0x138bca, 0x5e518, 0x7da4e2, 0x100001, 0x1300003, 0x400004, 0x2a0000a, 0x9800018, 0x200022, 0x27100061, 0x293000a3},
    {0x0, 0x0, 0x0, 0x7, 0xb, 0xf, 0x27, 0x4e, 0x9c, 0xd, 0x1e1, 0x6ed, 0x415, 0x82a, 0x1054, 0x4cc3, 0xddc7, 0x4e63, 0x3eb8b, 0x7d716, 0xfae2c, 0x8cee9, 0x9d875, 0x55bf89, 0x100111, 0x200222, 0x400444, 0x7f00fff, 0xab01aab, 0xdf02ddf, 0x21706117, 0x42e0c22e},
    {0x0, 0x1, 0x0, 0x6, 0x7, 0x7, 0x2d, 0x37, 0xb4, 0x6, 0x3d3, 0x4a3, 0x451, 0x1cf3, 0x1144, 0x5b6e, 0x28a7, 0xe797, 0x379fd, 0x18a07, 0xde7f4, 0x8dbe6, 0x23f223, 0x7b5253, 0x101101, 0x1303303, 0x404404, 0x6e0ee0e, 0x6716717, 0x5725727, 0x2bd6bd6d, 0x3c7bc7b7},
    {0x0, 0x0, 0x1, 0x2, 0xb, 0xc, 0x33, 0x66, 0xbf, 0x17e, 0x3d, 0x434, 0x505, 0xa0a, 0x111, 0x222, 0xc777, 0x5c9c, 0x28ebf, 0x51d7e, 0xcb443, 0x196886, 0x170dc9, 0x7f98e4, 0x110011, 0x220022, 0x1550055, 0x2aa00aa, 0xaab01ab, 0xeec02ec, 0x34430743, 0x68860e86},
    {0x0, 0x1, 0x2, 0x7, 0x8, 0x11, 0x2b, 0x3, 0x85, 0x10a, 0x391, 0x428, 0x8d3, 0x153d, 0x5, 0x400f, 0x801e, 0x1c033, 0x20078, 0x440f5, 0xac1c7, 0xc28f, 0x214791, 0x428f22, 0xe459d5, 0x10a3c88, 0x234bb9f, 0x54fa1c9, 0x11, 0x10000033, 0x20000066, 0x700000ff},
    {0x0, 0x1, 0x2, 0x2, 0x5, 0x19, 0x14, 0x9, 0x9b, 0x136, 0x15a, 0x23d, 0xd71, 0x8f4, 0x41, 0x40c3, 0x8186, 0x828a, 0x14555, 0x64e79, 0x51554, 0x262c9, 0x26a75b, 0x4d4eb6, 0x57d3da, 0x8dc57d, 0x35b0131, 0x23715f4, 0x1001, 0x10003003, 0x20006006, 0x2000a00a},
    {0x0, 0x0, 0x2, 0x3, 0xd, 0x1c, 0x2e, 0xf, 0x1e, 0x122, 0x1e9, 0x63b, 0xf54, 0x157a, 0x55, 0xaa, 0x81fe, 0xc257, 0x346f9, 0x70c0c, 0xb9cb6, 0x3e983, 0x7d306, 0x48750a, 0x7b769d, 0x18d9ba7, 0x3d34244, 0x55069b2, 0x1111, 0x2222, 0x20006666, 0x3000bbbb},
    {0x0, 0x1, 0x1, 0x1, 0x0, 0x6, 0x16, 0x11, 0xb3, 0xd5, 0x19, 0x110, 0x146, 0xe36, 0x101, 0x4303, 0x4505, 0x4909, 0x1010, 0x1a626, 0x5d656, 0x4d191, 0x2d72b3, 0x3797d5, 0x25d19, 0x4d1910, 0x40d746, 0x3a39836, 0x10001, 0x10030003, 0x10050005, 0x10090009},
    {0x0, 0x1, 0x0, 0x2, 0x3, 0x12, 0x3a, 0x1d, 0xa7, 0x74, 0x1d2, 0x77, 0xb4a, 0x18e2, 0x151, 0x43f3, 0x544, 0x882a, 0xd6e3, 0x4bd92, 0xee35a, 0x7f30d, 0x281517, 0x1fcc34, 0x707e72, 0x1725c7, 0x2ceb76a, 0x60dec42, 0x11101, 0x10033303, 0x44404, 0x200aaa0a},
    {0x0, 0x1, 0x0, 0x3, 0x9, 0xd, 0x3b, 0x27, 0xe9, 0x9c, 0xd1, 0x7ef, 0x3e3, 0x1381, 0x415, 0x4c3f, 0x1054, 0xec97, 0x225ed, 0x3f649, 0xf2a27, 0xb544b, 0x3dfcdd, 0x2d512c, 0x275e85, 0x1c4b2a3, 0xb62e5f, 0x42d8195, 0x100111, 0x10300333, 0x400444, 0x30b00bbb},
    {0x0, 0x0, 0x1, 0x4, 0xd, 0x4, 0x16, 0x2b, 0x56, 0x7, 0x3f4, 0x5ef, 0x7cc, 0x38a, 0x445, 0x88a, 0x5551, 0x1333c, 0x33329, 0x199b4, 0x4cc8e, 0x844c7, 0x10898e, 0x957db, 0xe33524, 0x14f3d93, 0x1a98bfc, 0x64e722, 0x101011, 0x202022, 0x10505055, 0x40c0c0cc},
    {0x0, 0x1, 0x1, 0x5, 0xd, 0x8, 0x1c, 0x39, 0xcb, 0x5d, 0x395, 0x405, 0x2e8, 0x2fc, 0x541, 0x4fc3, 0x5045, 0x17a4d, 0x32e5d, 0x28228, 0x63b5c, 0xc33f9, 0x34540b, 0x1cfc1d, 0xfd63d5, 0x13e5c45, 0xe7e0e8, 0x1e91fc, 0x111001, 0x10333003, 0x10555005, 0x50ddd00d},
    {0x0, 0x0, 0x3, 0x3, 0x0, 0x12, 0x0, 0x3f, 0x7e, 0x13d, 0x39, 0x3f0, 0xd6e, 0xfc0, 0x555, 0xaaa, 0xdaab, 0xe557, 0x5550, 0x4755a, 0x15540, 0xda573, 0x1b4ae6, 0x407a59, 0x1bc40d, 0xda5730, 0x335b3b6, 0x3695cc0, 0x111111, 0x222222, 0x30777777, 0x30bbbbbb},
    {0x0, 0x1, 0x2, 0x3, 0x5, 0x16, 0x38, 0x41, 0x43, 0x86, 0x34b, 0x7d5, 0x6b6, 0x278, 0x1001, 0x7003, 0xe006, 0x700b, 0x1015, 0x6e036, 0x98078, 0x1c50c1, 0x4f143, 0x9e286, 0xa6774b, 0x1081fd5, 0xc6f6b6, 0x2b9a278, 0x1000001, 0x13000003, 0x26000006, 0x3b00000b},
    {0x0, 0x0, 0x2, 0x1, 0x3, 0xd, 0xb, 0x4b, 0x96, 0xba, 0x293, 0x5ed, 0xcdf, 0x15c5, 0x1045, 0x208a, 0xe19e, 0xd26d, 0x1f49f, 0x19bd9, 0x663a7, 0x1e4127, 0x3c824e, 0x586d2, 0xcc481f, 0x1a6d119, 0x2fd6863, 0x6208391, 0x1001011, 0x2002022, 0x26006066, 0x19009099},
    {0x0, 0x1, 0x1, 0x2, 0xe, 0x0, 0x38, 0x53, 0x75, 0x19f, 0x33e, 0x142, 0xa60, 0x508, 0x1105, 0x730f, 0x1511, 0x2a22, 0x27e66, 0x220a0, 0x9f998, 0x19209f, 0xb61a1, 0x5da2e3, 0xbb45c6, 0xcdce4a, 0x32413e0, 0x3373928, 0x1010011, 0x13030033, 0x15050055, 0x2a0a00aa},
    {0x0, 0x1, 0x2, 0x6, 0x5, 0x11, 0x6, 0x55, 0x7f, 0xfe, 0x56, 0x6d1, 0x725, 0x17be, 0x1111, 0x7333, 0xe666, 0x16eee, 0x445, 0x76221, 0x5a226, 0x18c885, 0x9598f, 0x12b31e, 0x54f736, 0x15762c1, 0xad5075, 0x6a0925e, 0x1010101, 0x13030303, 0x26060606, 0x6e0e0e0e},
    {0x0, 0x1, 0x0, 0x4, 0x4, 0xe, 0x15, 0x65, 0x2f, 0x194, 0xbc, 0x5c4, 0x9d6, 0x1461, 0x1411, 0x7c33, 0x5044, 0x1f0cc, 0x1154, 0x1dace, 0x14105, 0x148ab5, 0x1d9fdf, 0x522ad4, 0x767f7c, 0x19a8184, 0x38e3cb6, 0x68ca671, 0x1100101, 0x13300303, 0x4400404, 0x4cc00c0c},
    {0x0, 0x0, 0x1, 0x2, 0xe, 0x1, 0x14, 0x6f, 0xde, 0x153, 0x2a6, 0x3ea, 0xd0f, 0x168c, 0x1455, 0x28aa, 0x501, 0xa02, 0x21e06, 0x2def5, 0x10144, 0x16ba63, 0x2d74c6, 0x6c53ef, 0xd8a7de, 0x69e862, 0x2e1f603, 0x61fd77c, 0x1101111, 0x2202222, 0x15505555, 0x2aa0aaaa},
    {0x0, 0x0, 0x3, 0x5, 0x0, 0x9, 0x20, 0x71, 0xe2, 0xd7, 0xbd, 0x710, 0x959, 0x260, 0x1501, 0x2a02, 0xab07, 0x1a90d, 0x15010, 0x5d29, 0xfe060, 0x1065f1, 0x20cbe2, 0x1139d7, 0x72ddbd, 0x1065f10, 0x3bff459, 0x215c260, 0x1110001, 0x2220002, 0x37770007, 0x5ddd000d},
    {0x0, 0x0, 0x3, 0x1, 0x5, 0xd, 0x5, 0x77, 0xee, 0xc5, 0x34f, 0x45b, 0xa73, 0x1eeb, 0x1515, 0x2a2a, 0xab6b, 0xfdbd, 0x5011, 0x10b49, 0x44401, 0x11975b, 0x232eb6, 0x14e481, 0xbd2d83, 0x1eebf87, 0x3499b8f, 0x4921cf7, 0x1110111, 0x2220222, 0x37770777, 0x19990999},
    {0x0, 0x1, 0x3, 0x3, 0xb, 0x12, 0x3, 0x7d, 0x7, 0xf3, 0x2ef, 0x13f, 0x18a, 0x1e47, 0x1551, 0x7ff3, 0xaab7, 0x557b, 0x3006b, 0x75592, 0x5abb3, 0x13a96d, 0x14fbb7, 0x1a5e03, 0xc9b0df, 0xf3260f, 0x128e9aa, 0x4bea0f7, 0x1111101, 0x13333303, 0x37777707, 0x3bbbbb0b},
    {0x0, 0x1, 0x0, 0x6, 0x8, 0xa, 0x7, 0x22, 0x1d, 0x127, 0x74, 0x6a6, 0x938, 0x972, 0x13, 0x2f1a, 0x151, 0x103f3, 0x544, 0x60d6e, 0x81f98, 0xa220a, 0x752f7, 0x228002, 0x1c4a8d, 0x124df97, 0x712a34, 0x6abeb46, 0x926fcb8, 0x95390d2, 0x4756e3, 0x2f94823a},
    {0x0, 0x1, 0x1, 0x7, 0x8, 0x1f, 0x2e, 0x6, 0x2b, 0x17d, 0x187, 0x689, 0xbe8, 0x1959, 0x2002, 0x137a, 0x445, 0x10ccf, 0x11551, 0x73ff3, 0x86678, 0x1ff303, 0x2fa256, 0x43b1e, 0x2fe347, 0x17025c9, 0x1906e5b, 0x6b0b2ed, 0xb812e48, 0x19b2ee7d, 0x219be88a, 0x1111e812},
    {0x0, 0x1, 0x0, 0x5, 0xe, 0x1a, 0x12, 0x15, 0x2d, 0x177, 0xb4, 0x4f1, 0xd56, 0x1c42, 0x1bca, 0x1c9, 0x451, 0x10cf3, 0x1144, 0x5379d, 0xe7e7e, 0x1ae51a, 0x1359f2, 0x177885, 0x29ecbd, 0x17a35c7, 0xa7b2f4, 0x4c13ba1, 0xd25c5b6, 0x1cbfe0e2, 0x1ab63dea, 0x3e6cb19},
    {0x0, 0x1, 0x0, 0x5, 0xf, 0x1a, 0x24, 0x1f, 0x4d, 0x1d7, 0x134, 0x611, 0x85b, 0x1582, 0x3fd4, 0x3edb, 0x1051, 0x130f3, 0x4144, 0x5d39d, 0xef62f, 0x19ad1a, 0x225f24, 0x16deaf, 0x5996dd, 0x1eabb67, 0x1665b74, 0x6f37b41, 0x9d93b2b, 0x16d42d22, 0x38313794, 0x351255ab},
    {0x0, 0x0, 0x3, 0x4, 0x6, 0x11, 0x6, 0x20, 0x5f, 0xbe, 0x29d, 0x784, 0x232, 0x1f4f, 0x1002, 0x460, 0x1155, 0x22aa, 0x376ab, 0x4cffc, 0x772ae, 0x122ea5, 0x232be, 0x2a8020, 0x4be913, 0x97d226, 0x2f39f79, 0x770ecd4, 0x306e75a, 0x1c885a43, 0x154232aa, 0xc89abe0},
    {0x0, 0x1, 0x2, 0x6, 0x2, 0x19, 0x1, 0x69, 0x63, 0x1a5, 0x34a, 0x452, 0x4f6, 0x102b, 0x19a3, 0x4f5b, 0x1405, 0x13c0f, 0x2781e, 0x6d836, 0x3685a, 0x1a74dd, 0x41545, 0x64374d, 0x70b8ef, 0x191c931, 0x3239262, 0x4a6551a, 0x5eaff2e, 0x13e9ec87, 0x1d5e832f, 0x40902e37},
    {0x0, 0x1, 0x0, 0x0, 0x9, 0x1d, 0x3d, 0x26, 0x65, 0x1af, 0x194, 0x328, 0xc1d, 0x1529, 0x2c69, 0x197e, 0x1411, 0x13c33, 0x5044, 0xa088, 0x8f589, 0x1e27ed, 0x3b23ad, 0x2ef2c6, 0x76d335, 0x19b755f, 0x1db4cd4, 0x3b699a8, 0xdad79cd, 0x16ac53b9, 0x2b189ef9, 0x1285169e},
    {0x0, 0x0, 0x1, 0x1, 0x8, 0xb, 0x14, 0x3b, 0x69, 0xd2, 0xcd, 0x221, 0xdd8, 0x5d3, 0x974, 0x7c3, 0x1441, 0x2882, 0x14545, 0x1b649, 0x9e618, 0x916eb, 0x100554, 0x32727b, 0x7a6f29, 0xf4de52, 0x93d38d, 0x2a91661, 0xc758bd8, 0x7102d13, 0xcd48474, 0xd814b03},
    {0x0, 0x1, 0x2, 0x7, 0x0, 0x1, 0x8, 0x18, 0x71, 0x193, 0x326, 0x5df, 0x710, 0xf51, 0x17c8, 0x2418, 0x1501, 0x13f03, 0x27e06, 0x7c30f, 0x15010, 0x3b521, 0xde848, 0x137898, 0x62a471, 0x1a7ec93, 0x34fd926, 0x5385edf, 0x62a4710, 0xd362a51, 0x13bc3fc8, 0x2c6d5c18},
    {0x0, 0x1, 0x3, 0x7, 0xb, 0x7, 0x32, 0x56, 0x87, 0x89, 0x95, 0xad, 0x6c1, 0x1475, 0xa5e, 0x3f22, 0x4015, 0x1c03f, 0x2c06b, 0x4c0c3, 0xdc1c7, 0xec2cb, 0x2e869a, 0x638eee, 0xe6dfeb, 0x2b603d, 0x1b01f91, 0x286e0c9, 0x37061d5, 0x196be2f1, 0x18cc46c6, 0x10e3318a},
    {0x0, 0x1, 0x2, 0x7, 0xa, 0x14, 0x2a, 0x7f, 0x8d, 0x97, 0x12e, 0xcb, 0x7a2, 0xf44, 0x1d92, 0x47b, 0x4051, 0x1c0f3, 0x381e6, 0x4c33f, 0xc873a, 0x190e74, 0x309c4a, 0x40f0cf, 0xee3a1d, 0x324e27, 0x649c4e, 0x2fb76bb, 0x24e0502, 0x49c0a04, 0xae46032, 0x254c1d0b},
    {0x0, 0x0, 0x2, 0x0, 0xa, 0x16, 0x3, 0x5, 0xa9, 0x152, 0x1f6, 0x548, 0x48a, 0xa46, 0x28bb, 0x538d, 0x4441, 0x8882, 0x39986, 0x22208, 0xcee9a, 0x1b55b6, 0x12dc83, 0x2675c5, 0xc58fe9, 0x18b1fd2, 0x9d2076, 0x62c7f48, 0x1ff9e0a, 0x7423c6, 0x332d6a7b, 0x641444cd},
    {0x0, 0x0, 0x0, 0x4, 0x9, 0x1c, 0xa, 0x9, 0xc3, 0x186, 0x30c, 0x114, 0x3eb, 0xd44, 0x3d5e, 0x6e5b, 0x5005, 0xa00a, 0x14014, 0x7c03c, 0xed07d, 0x10c0cc, 0x1c2162, 0x23d2ad, 0xaff6cf, 0x15fed9e, 0x2bfdb3c, 0x3c06d44, 0x62f2c47, 0x4c1d854, 0x25dde826, 0x5b2b2737},
    {0x0, 0x1, 0x2, 0x5, 0x1, 0x3, 0x12, 0x72, 0xcf, 0x51, 0xa2, 0x8b, 0xd3f, 0x1bb1, 0x2cae, 0x32ce, 0x5055, 0x1f0ff, 0x3e1fe, 0x693a9, 0x45505, 0x9fa5f, 0x3b0ba, 0x41909a, 0xa059c3, 0xe0ea45, 0x1c1d48a, 0x223f0d7, 0xba5c5f3, 0x16ebd225, 0x31535f76, 0x15748696},
    {0x0, 0x1, 0x0, 0x3, 0x9, 0x14, 0x39, 0x3d, 0xe7, 0x29, 0x39c, 0x511, 0xaf, 0x50c, 0x158f, 0x5853, 0x5415, 0x1fc3f, 0x15054, 0x15c97, 0xeb5ed, 0x1a93a4, 0x21320d, 0x1b6d99, 0x87340b, 0x895c1d, 0x21cd02c, 0x6b0fc45, 0x5cdd4e3, 0xe8911fc, 0x4e65743, 0x60ad802f},
    {0x0, 0x0, 0x2, 0x2, 0x8, 0x1f, 0x1c, 0x3e, 0xf5, 0x1ea, 0x3e, 0x442, 0xf8, 0xb93, 0x2a6c, 0x50e6, 0x5511, 0xaa22, 0x3fe66, 0x2aa, 0xff998, 0x13f0cf, 0xfe99c, 0x182d5e, 0x912fa5, 0x1225f4a, 0x166e1de, 0x7ab2262, 0x59b8778, 0x349bda3, 0x3794d0ac, 0x6a4f4086},
    {0x0, 0x1, 0x1, 0x3, 0x7, 0x3, 0x3c, 0x30, 0xff, 0x11, 0x233, 0x255, 0x6bb, 0xf67, 0x413, 0x7fbc, 0x6bb0, 0x1e00f, 0x101, 0x40303, 0x40505, 0xc0b0b, 0x1c1717, 0xc2323, 0xf07c7c, 0xc0b0b0, 0x3fdfeff, 0x461311, 0x8ca3533, 0x95e5f55, 0x1afaadbb, 0x3db34867},
    {0x0, 0x1, 0x3, 0x6, 0xe, 0x0, 0x12, 0x1e, 0x55, 0x1b, 0x22d, 0x641, 0xc82, 0x1d32, 0x360, 0x2346, 0x30b2, 0xb607, 0x145, 0x403cf, 0xc06db, 0x180db6, 0x3819e6, 0x28a0, 0x48479a, 0x78bb66, 0x1550441, 0x6e97b7, 0x8b3b8d9, 0x1909e605, 0x3213cc0a, 0x74fab77a},
    {0x0, 0x0, 0x0, 0x2, 0x1, 0x1a, 0x2f, 0x3a, 0xc1, 0x21, 0x42, 0x84, 0x54a, 0x31, 0x337a, 0x538f, 0x63fa, 0x1bbe1, 0x401, 0x802, 0x1004, 0x8280a, 0x44411, 0x68e83a, 0xbdbc6f, 0xeae8ba, 0x30305c1, 0x8c8621, 0x1190c42, 0x2321884, 0x157d3d4a, 0x44e431},
    {0x0, 0x1, 0x0, 0x3, 0x2, 0x1d, 0x20, 0x5a, 0xc5, 0x2d, 0x277, 0xb4, 0x71f, 0x68a, 0x3c81, 0x4ee0, 0xaa22, 0x1ba59, 0x451, 0x40cf3, 0x1144, 0xc2e7b, 0x84db2, 0x74f8ad, 0x819e60, 0x16b53fa, 0x31378d5, 0xbc1fbd, 0x9c420c7, 0x2f07ef4, 0x1c24dd2f, 0x1ab9c4aa},
    {0x0, 0x0, 0x0, 0x5, 0xf, 0x2, 0x36, 0x19, 0xff, 0x33, 0x66, 0xcc, 0xb67, 0x1c31, 0x206, 0x653a, 0x291b, 0x1dc11, 0x505, 0xa0a, 0x1414, 0x143939, 0x3c6363, 0x8aaaa, 0xd9afae, 0x66fffd, 0x3fa0503, 0xc6f5ff, 0x18debfe, 0x31bd7fc, 0x2dea8dfb, 0x700839f5},
    {0x0, 0x0, 0x1, 0x4, 0x2, 0xb, 0x3, 0x6e, 0xcd, 0x59, 0xb2, 0x33d, 0xbac, 0x122, 0x1f03, 0x10ab, 0xeefe, 0x1fa35, 0x1141, 0x2282, 0x45445, 0x10cf0c, 0x93692, 0x2e91eb, 0x86383, 0x1b6356e, 0x3296f8d, 0x1435d19, 0x286ba32, 0xc4e297d, 0x2f179cac, 0x6b36ba2},
    {0x0, 0x0, 0x1, 0x6, 0xf, 0xa, 0xc, 0x66, 0x95, 0x5f, 0xbe, 0x323, 0xf3a, 0x1895, 0x1da6, 0xc44, 0xfe62, 0x15e53, 0x1155, 0x22aa, 0x45401, 0x18ed56, 0x3de953, 0x2a82a2, 0x349abc, 0x196b29e, 0x24cbed1, 0x15b1613, 0x2b62c26, 0xc374e5f, 0x3d02c4f2, 0x61e8b3d1},
    {0x0, 0x1, 0x2, 0x5, 0x9, 0x1a, 0x3d, 0x18, 0x5b, 0x69, 0x2bb, 0x576, 0x885, 0x17b1, 0x3c2a, 0x6975, 0x158, 0xc023, 0x1441, 0x43cc3, 0x87986, 0x14e74d, 0x25f259, 0x6b46ba, 0xf23b3d, 0x6bc698, 0x17c8b9b, 0x18bac29, 0xa9cf47b, 0x1539e8f6, 0x23f87dc5, 0x5d6c0ff1},
    {0x0, 0x0, 0x0, 0x4, 0x9, 0x11, 0x27, 0x40, 0xdd, 0x6f, 0xde, 0x1bc, 0xac4, 0x17e7, 0x297f, 0x592d, 0xac40, 0x1fd1b, 0x1455, 0x28aa, 0x5154, 0x10f3fc, 0x25f3ad, 0x47dba5, 0x9bf24b, 0x10f3fc0, 0x36ec839, 0x193fae3, 0x327f5c6, 0x64feb8c, 0x2ad03c94, 0x5c3383cb},
    {0x0, 0x0, 0x1, 0x3, 0xd, 0x13, 0x29, 0x16, 0x9, 0x77, 0xee, 0x3ab, 0x521, 0x1f63, 0x2f09, 0x42ef, 0x11c2, 0x66cf, 0x1515, 0x2a2a, 0x44141, 0xc9797, 0x35b8b9, 0x4fcccf, 0xa35a5d, 0x53a5ae, 0x31a8bd, 0x1f0f7db, 0x3e1efb6, 0xe3328b7, 0x1596a6b5, 0x7ebbebdf},
    {0x0, 0x0, 0x2, 0x3, 0x3, 0x1f, 0x2e, 0x58, 0x1, 0x7d, 0xfa, 0x50e, 0x56f, 0x157, 0x34eb, 0x4e06, 0x95f8, 0x7f7d, 0x1551, 0x2aa2, 0x87fe6, 0xc957b, 0xd6ae3, 0x7f3f0f, 0xbf2b0e, 0x16e0358, 0x114451, 0x1d8e3ed, 0x3b1c7da, 0x14d2486e, 0x14ae3b5f, 0x7e71ae7},
};

static constexpr int FAURE03_XOR_SIZE = 21;
static constexpr uint32_t faure03_xors[3][FAURE03_XOR_SIZE] = {
    {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},
    {0x0, 0x1, 0x8, 0x1, 0x20, 0xe5, 0x38, 0x39d, 0x19a0, 0x1, 0x4ce8, 0x26725, 0x4d00, 0x99cee, 0x44c888, 0x10d49d, 0x115df00, 0x7b25f51, 0x99c8, 0x1719ab85, 0xb8c45970},
    {0x0, 0x2, 0x5, 0x2, 0x3e, 0x9e, 0x1d, 0x626, 0x1004, 0x2, 0x99ce, 0x18086, 0x99fe, 0x129fda, 0x2f769e, 0x8bb86, 0x1d8d2ee, 0x4cfa6d2, 0x4ce5, 0x2e305626, 0x737aa4b4},
};

static constexpr int FAURE05_XOR_SIZE = 14;
static constexpr uint32_t faure05_xors[5][FAURE05_XOR_SIZE] = {
    {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},
    {0x0, 0x1, 0xe, 0x56, 0x270, 0x1, 0xc3e, 0xab0f, 0x41a90, 0x1dc43d, 0x186e, 0x95a1bf, 0x828e1dd, 0x3220140e},
    {0x0, 0x2, 0x15, 0x2b, 0x19e, 0x2, 0x1875, 0x10098, 0x20e44, 0x13c4de, 0x30d5, 0x12b1288, 0xc3f039c, 0x1921652b},
    {0x0, 0x3, 0x6, 0x75, 0x126, 0x3, 0x24af, 0x4995, 0x59604, 0xe0e5f, 0xc36, 0x1bf5e65, 0x380682c, 0x4422222f},
    {0x0, 0x4, 0x13, 0x40, 0xea, 0x4, 0x30ec, 0xe866, 0x30f74, 0xb3320, 0x24a3, 0x254e7b6, 0xb132e21, 0x25587f8e},
};

static constexpr int FAURE07_XOR_SIZE = 12;
static constexpr uint32_t faure07_xors[7][FAURE07_XOR_SIZE] = {
    {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},
    {0x0, 0x1, 0x14, 0xb0, 0x5af, 0x34e3, 0x1cb90, 0x1, 0xc9104, 0xfb53a1, 0x8a3ac06, 0x476bfa7b},
    {0x0, 0x2, 0x1f, 0x135, 0x23c, 0x1f84, 0x16d80, 0x2, 0x1921ff, 0x1858e76, 0xf2afd4d, 0x1c1402ad},
    {0x0, 0x3, 0x2f, 0x6f, 0x72d, 0xb88, 0x1377f, 0x3, 0x25b2ff, 0x24e9e07, 0x572df6e, 0x5a2c60f1},
    {0x0, 0x4, 0xc, 0xfd, 0x32f, 0x3a39, 0xfca0, 0x4, 0x3243fd, 0x96cc6e, 0xc6b49c0, 0x2801b132},
    {0x0, 0x5, 0x18, 0x45, 0x8e8, 0x2915, 0xa1df, 0x5, 0x3ed4f9, 0x12d984d, 0x36319c9, 0x6feb4dce},
    {0x0, 0x6, 0x29, 0xe6, 0x451, 0x18f6, 0x7c77, 0x6, 0x4b65fa, 0x20338de, 0xb4a46da, 0x363df6be},
};

static constexpr int FAURE011_XOR_SIZE = 10;
static constexpr uint32_t faure011_xors[11][FAURE011_XOR_SIZE] = {
    {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},
    {0x0, 0x1, 0x20, 0x1c4, 0x175f, 0x12825, 0x10815f, 0xc0544e, 0x9d60449, 0x7d15b085},
    {0x0, 0x2, 0x33, 0x34c, 0x2e4d, 0x25283, 0x3ec9e, 0x5c10b0, 0x6d65a04, 0x653d9e07},
    {0x0, 0x3, 0x44, 0x488, 0x616, 0x1167c, 0x130598, 0x126bfde, 0x27b63bb, 0x471683f6},
    {0x0, 0x4, 0x5e, 0xcf, 0x1b93, 0x21df8, 0x57667, 0xb029c9, 0xc06fdb9, 0x2f60b5a5},
    {0x0, 0x5, 0x76, 0x1fe, 0x30e8, 0xbeea, 0x15f86b, 0x3f5317, 0x89939a8, 0xfd5a44d},
    {0x0, 0x6, 0x13, 0x36c, 0xc8e, 0x1db7d, 0x9c8d8, 0xfb6a05, 0x524bac8, 0x82b88868},
    {0x0, 0x7, 0x27, 0x509, 0x20fd, 0x90f9, 0x16b9b2, 0x9471cd, 0x181d0cb, 0x6ebc105e},
    {0x0, 0x8, 0x39, 0x13a, 0x3494, 0x1c19f, 0xbd507, 0x316f6c, 0xa9f7beb, 0x549ad223},
    {0x0, 0x9, 0x54, 0x2ce, 0x13eb, 0x521e, 0x1a0569, 0xe376ec, 0x7f33bd8, 0x3fdc17e2},
    {0x0, 0xa, 0x6d, 0x42a, 0x2719, 0x15d1e, 0xde7d4, 0x70c7f2, 0x41d9001, 0x2354d3fa},
};

}  // namespace sampling

#endif  // SAMPLING_XOR_VALUES_H