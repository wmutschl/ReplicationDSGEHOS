function nXI3min = Xi_Gaussian_prodmom3_orderApp3_nx3_nu1(arg)
SIGe_1 = arg(1);
E_XF1_1 = arg(2);
E_XF1_2 = arg(3);
E_XF1_3 = arg(4);
E_XF2_1 = arg(5);
E_XF2_2 = arg(6);
E_XF2_3 = arg(7);
E_XF2_4 = arg(8);
E_XF2_5 = arg(9);
E_XF2_6 = arg(10);
E_XF3_1 = arg(11);
E_XF3_2 = arg(12);
E_XF3_3 = arg(13);
E_XF3_4 = arg(14);
E_XF3_5 = arg(15);
E_XF3_6 = arg(16);
E_XF3_7 = arg(17);
E_XF3_8 = arg(18);
E_XF3_9 = arg(19);
E_XF3_10 = arg(20);
E_XF4_1 = arg(21);
E_XF4_2 = arg(22);
E_XF4_3 = arg(23);
E_XF4_4 = arg(24);
E_XF4_5 = arg(25);
E_XF4_6 = arg(26);
E_XF4_7 = arg(27);
E_XF4_8 = arg(28);
E_XF4_9 = arg(29);
E_XF4_10 = arg(30);
E_XF4_11 = arg(31);
E_XF4_12 = arg(32);
E_XF4_13 = arg(33);
E_XF4_14 = arg(34);
E_XF4_15 = arg(35);
E_XF5_1 = arg(36);
E_XF5_2 = arg(37);
E_XF5_3 = arg(38);
E_XF5_4 = arg(39);
E_XF5_5 = arg(40);
E_XF5_6 = arg(41);
E_XF5_7 = arg(42);
E_XF5_8 = arg(43);
E_XF5_9 = arg(44);
E_XF5_10 = arg(45);
E_XF5_11 = arg(46);
E_XF5_12 = arg(47);
E_XF5_13 = arg(48);
E_XF5_14 = arg(49);
E_XF5_15 = arg(50);
E_XF5_16 = arg(51);
E_XF5_17 = arg(52);
E_XF5_18 = arg(53);
E_XF5_19 = arg(54);
E_XF5_20 = arg(55);
E_XF5_21 = arg(56);
E_XF6_1 = arg(57);
E_XF6_2 = arg(58);
E_XF6_3 = arg(59);
E_XF6_4 = arg(60);
E_XF6_5 = arg(61);
E_XF6_6 = arg(62);
E_XF6_7 = arg(63);
E_XF6_8 = arg(64);
E_XF6_9 = arg(65);
E_XF6_10 = arg(66);
E_XF6_11 = arg(67);
E_XF6_12 = arg(68);
E_XF6_13 = arg(69);
E_XF6_14 = arg(70);
E_XF6_15 = arg(71);
E_XF6_16 = arg(72);
E_XF6_17 = arg(73);
E_XF6_18 = arg(74);
E_XF6_19 = arg(75);
E_XF6_20 = arg(76);
E_XF6_21 = arg(77);
E_XF6_22 = arg(78);
E_XF6_23 = arg(79);
E_XF6_24 = arg(80);
E_XF6_25 = arg(81);
E_XF6_26 = arg(82);
E_XF6_27 = arg(83);
E_XF6_28 = arg(84);
E_XS1_1 = arg(85);
E_XS1_2 = arg(86);
E_XS1_3 = arg(87);
E_XS2_1 = arg(88);
E_XS2_2 = arg(89);
E_XS2_3 = arg(90);
E_XS2_4 = arg(91);
E_XS2_5 = arg(92);
E_XS2_6 = arg(93);
E_XS3_1 = arg(94);
E_XS3_2 = arg(95);
E_XS3_3 = arg(96);
E_XS3_4 = arg(97);
E_XS3_5 = arg(98);
E_XS3_6 = arg(99);
E_XS3_7 = arg(100);
E_XS3_8 = arg(101);
E_XS3_9 = arg(102);
E_XS3_10 = arg(103);
E_XF1_XS1_1 = arg(104);
E_XF1_XS1_2 = arg(105);
E_XF1_XS1_3 = arg(106);
E_XF1_XS1_4 = arg(107);
E_XF1_XS1_5 = arg(108);
E_XF1_XS1_6 = arg(109);
E_XF1_XS1_7 = arg(110);
E_XF1_XS1_8 = arg(111);
E_XF1_XS1_9 = arg(112);
E_XF2_XS1_1 = arg(113);
E_XF2_XS1_2 = arg(114);
E_XF2_XS1_3 = arg(115);
E_XF2_XS1_4 = arg(116);
E_XF2_XS1_5 = arg(117);
E_XF2_XS1_6 = arg(118);
E_XF2_XS1_7 = arg(119);
E_XF2_XS1_8 = arg(120);
E_XF2_XS1_9 = arg(121);
E_XF2_XS1_10 = arg(122);
E_XF2_XS1_11 = arg(123);
E_XF2_XS1_12 = arg(124);
E_XF2_XS1_13 = arg(125);
E_XF2_XS1_14 = arg(126);
E_XF2_XS1_15 = arg(127);
E_XF2_XS1_16 = arg(128);
E_XF2_XS1_17 = arg(129);
E_XF2_XS1_18 = arg(130);
E_XF1_XS2_1 = arg(131);
E_XF1_XS2_2 = arg(132);
E_XF1_XS2_3 = arg(133);
E_XF1_XS2_4 = arg(134);
E_XF1_XS2_5 = arg(135);
E_XF1_XS2_6 = arg(136);
E_XF1_XS2_7 = arg(137);
E_XF1_XS2_8 = arg(138);
E_XF1_XS2_9 = arg(139);
E_XF1_XS2_10 = arg(140);
E_XF1_XS2_11 = arg(141);
E_XF1_XS2_12 = arg(142);
E_XF1_XS2_13 = arg(143);
E_XF1_XS2_14 = arg(144);
E_XF1_XS2_15 = arg(145);
E_XF1_XS2_16 = arg(146);
E_XF1_XS2_17 = arg(147);
E_XF1_XS2_18 = arg(148);
E_XF3_XS1_1 = arg(149);
E_XF3_XS1_2 = arg(150);
E_XF3_XS1_3 = arg(151);
E_XF3_XS1_4 = arg(152);
E_XF3_XS1_5 = arg(153);
E_XF3_XS1_6 = arg(154);
E_XF3_XS1_7 = arg(155);
E_XF3_XS1_8 = arg(156);
E_XF3_XS1_9 = arg(157);
E_XF3_XS1_10 = arg(158);
E_XF3_XS1_11 = arg(159);
E_XF3_XS1_12 = arg(160);
E_XF3_XS1_13 = arg(161);
E_XF3_XS1_14 = arg(162);
E_XF3_XS1_15 = arg(163);
E_XF3_XS1_16 = arg(164);
E_XF3_XS1_17 = arg(165);
E_XF3_XS1_18 = arg(166);
E_XF3_XS1_19 = arg(167);
E_XF3_XS1_20 = arg(168);
E_XF3_XS1_21 = arg(169);
E_XF3_XS1_22 = arg(170);
E_XF3_XS1_23 = arg(171);
E_XF3_XS1_24 = arg(172);
E_XF3_XS1_25 = arg(173);
E_XF3_XS1_26 = arg(174);
E_XF3_XS1_27 = arg(175);
E_XF3_XS1_28 = arg(176);
E_XF3_XS1_29 = arg(177);
E_XF3_XS1_30 = arg(178);
E_XF2_XS2_1 = arg(179);
E_XF2_XS2_2 = arg(180);
E_XF2_XS2_3 = arg(181);
E_XF2_XS2_4 = arg(182);
E_XF2_XS2_5 = arg(183);
E_XF2_XS2_6 = arg(184);
E_XF2_XS2_7 = arg(185);
E_XF2_XS2_8 = arg(186);
E_XF2_XS2_9 = arg(187);
E_XF2_XS2_10 = arg(188);
E_XF2_XS2_11 = arg(189);
E_XF2_XS2_12 = arg(190);
E_XF2_XS2_13 = arg(191);
E_XF2_XS2_14 = arg(192);
E_XF2_XS2_15 = arg(193);
E_XF2_XS2_16 = arg(194);
E_XF2_XS2_17 = arg(195);
E_XF2_XS2_18 = arg(196);
E_XF2_XS2_19 = arg(197);
E_XF2_XS2_20 = arg(198);
E_XF2_XS2_21 = arg(199);
E_XF2_XS2_22 = arg(200);
E_XF2_XS2_23 = arg(201);
E_XF2_XS2_24 = arg(202);
E_XF2_XS2_25 = arg(203);
E_XF2_XS2_26 = arg(204);
E_XF2_XS2_27 = arg(205);
E_XF2_XS2_28 = arg(206);
E_XF2_XS2_29 = arg(207);
E_XF2_XS2_30 = arg(208);
E_XF2_XS2_31 = arg(209);
E_XF2_XS2_32 = arg(210);
E_XF2_XS2_33 = arg(211);
E_XF2_XS2_34 = arg(212);
E_XF2_XS2_35 = arg(213);
E_XF2_XS2_36 = arg(214);
E_XF4_XS1_1 = arg(215);
E_XF4_XS1_2 = arg(216);
E_XF4_XS1_3 = arg(217);
E_XF4_XS1_4 = arg(218);
E_XF4_XS1_5 = arg(219);
E_XF4_XS1_6 = arg(220);
E_XF4_XS1_7 = arg(221);
E_XF4_XS1_8 = arg(222);
E_XF4_XS1_9 = arg(223);
E_XF4_XS1_10 = arg(224);
E_XF4_XS1_11 = arg(225);
E_XF4_XS1_12 = arg(226);
E_XF4_XS1_13 = arg(227);
E_XF4_XS1_14 = arg(228);
E_XF4_XS1_15 = arg(229);
E_XF4_XS1_16 = arg(230);
E_XF4_XS1_17 = arg(231);
E_XF4_XS1_18 = arg(232);
E_XF4_XS1_19 = arg(233);
E_XF4_XS1_20 = arg(234);
E_XF4_XS1_21 = arg(235);
E_XF4_XS1_22 = arg(236);
E_XF4_XS1_23 = arg(237);
E_XF4_XS1_24 = arg(238);
E_XF4_XS1_25 = arg(239);
E_XF4_XS1_26 = arg(240);
E_XF4_XS1_27 = arg(241);
E_XF4_XS1_28 = arg(242);
E_XF4_XS1_29 = arg(243);
E_XF4_XS1_30 = arg(244);
E_XF4_XS1_31 = arg(245);
E_XF4_XS1_32 = arg(246);
E_XF4_XS1_33 = arg(247);
E_XF4_XS1_34 = arg(248);
E_XF4_XS1_35 = arg(249);
E_XF4_XS1_36 = arg(250);
E_XF4_XS1_37 = arg(251);
E_XF4_XS1_38 = arg(252);
E_XF4_XS1_39 = arg(253);
E_XF4_XS1_40 = arg(254);
E_XF4_XS1_41 = arg(255);
E_XF4_XS1_42 = arg(256);
E_XF4_XS1_43 = arg(257);
E_XF4_XS1_44 = arg(258);
E_XF4_XS1_45 = arg(259);
nXI3min=zeros(1140,1);
nXI3min(2,1) = 2*SIGe_1^2;
nXI3min(15,1) = 3*E_XF1_1*SIGe_1^2;
nXI3min(16,1) = 3*E_XF1_2*SIGe_1^2;
nXI3min(17,1) = 3*E_XF1_3*SIGe_1^2;
nXI3min(20,1) = 2*E_XF1_1*SIGe_1^2;
nXI3min(21,1) = 2*E_XF1_2*SIGe_1^2;
nXI3min(22,1) = 2*E_XF1_3*SIGe_1^2;
nXI3min(23,1) = 2*E_XS1_1*SIGe_1^2;
nXI3min(24,1) = 2*E_XS1_2*SIGe_1^2;
nXI3min(25,1) = 2*E_XS1_3*SIGe_1^2;
nXI3min(26,1) = 2*E_XF2_1*SIGe_1^2;
nXI3min(27,1) = 2*E_XF2_2*SIGe_1^2;
nXI3min(28,1) = 2*E_XF2_3*SIGe_1^2;
nXI3min(29,1) = 2*E_XF2_4*SIGe_1^2;
nXI3min(30,1) = 2*E_XF2_5*SIGe_1^2;
nXI3min(31,1) = 2*E_XF2_6*SIGe_1^2;
nXI3min(35,1) = 12*SIGe_1^3;
nXI3min(48,1) = 3*E_XF2_1*SIGe_1^2;
nXI3min(49,1) = 3*E_XF2_2*SIGe_1^2;
nXI3min(50,1) = 3*E_XF2_3*SIGe_1^2;
nXI3min(63,1) = 3*E_XF2_2*SIGe_1^2;
nXI3min(64,1) = 3*E_XF2_4*SIGe_1^2;
nXI3min(65,1) = 3*E_XF2_5*SIGe_1^2;
nXI3min(77,1) = 3*E_XF2_3*SIGe_1^2;
nXI3min(78,1) = 3*E_XF2_5*SIGe_1^2;
nXI3min(79,1) = 3*E_XF2_6*SIGe_1^2;
nXI3min(90,1) = 3*E_XF1_XS1_1*SIGe_1^2;
nXI3min(91,1) = 3*E_XF1_XS1_4*SIGe_1^2;
nXI3min(92,1) = 3*E_XF1_XS1_7*SIGe_1^2;
nXI3min(102,1) = 3*E_XF1_XS1_2*SIGe_1^2;
nXI3min(103,1) = 3*E_XF1_XS1_5*SIGe_1^2;
nXI3min(104,1) = 3*E_XF1_XS1_8*SIGe_1^2;
nXI3min(113,1) = 3*E_XF1_XS1_3*SIGe_1^2;
nXI3min(114,1) = 3*E_XF1_XS1_6*SIGe_1^2;
nXI3min(115,1) = 3*E_XF1_XS1_9*SIGe_1^2;
nXI3min(123,1) = 3*E_XF3_1*SIGe_1^2;
nXI3min(124,1) = 3*E_XF3_2*SIGe_1^2;
nXI3min(125,1) = 3*E_XF3_3*SIGe_1^2;
nXI3min(132,1) = 3*E_XF3_2*SIGe_1^2;
nXI3min(133,1) = 3*E_XF3_4*SIGe_1^2;
nXI3min(134,1) = 3*E_XF3_5*SIGe_1^2;
nXI3min(140,1) = 3*E_XF3_3*SIGe_1^2;
nXI3min(141,1) = 3*E_XF3_5*SIGe_1^2;
nXI3min(142,1) = 3*E_XF3_6*SIGe_1^2;
nXI3min(147,1) = 3*E_XF3_4*SIGe_1^2;
nXI3min(148,1) = 3*E_XF3_7*SIGe_1^2;
nXI3min(149,1) = 3*E_XF3_8*SIGe_1^2;
nXI3min(153,1) = 3*E_XF3_5*SIGe_1^2;
nXI3min(154,1) = 3*E_XF3_8*SIGe_1^2;
nXI3min(155,1) = 3*E_XF3_9*SIGe_1^2;
nXI3min(158,1) = 3*E_XF3_6*SIGe_1^2;
nXI3min(159,1) = 3*E_XF3_9*SIGe_1^2;
nXI3min(160,1) = 3*E_XF3_10*SIGe_1^2;
nXI3min(165,1) = 15*E_XF1_1*SIGe_1^3;
nXI3min(168,1) = 15*E_XF1_2*SIGe_1^3;
nXI3min(170,1) = 15*E_XF1_3*SIGe_1^3;
nXI3min(172,1) = 8*SIGe_1^3;
nXI3min(185,1) = 10*E_XF1_1*SIGe_1^3;
nXI3min(186,1) = 10*E_XF1_2*SIGe_1^3;
nXI3min(187,1) = 10*E_XF1_3*SIGe_1^3;
nXI3min(189,1) = 2*E_XF2_1*SIGe_1^2;
nXI3min(190,1) = 2*E_XF2_2*SIGe_1^2;
nXI3min(191,1) = 2*E_XF2_3*SIGe_1^2;
nXI3min(192,1) = 2*E_XF1_XS1_1*SIGe_1^2;
nXI3min(193,1) = 2*E_XF1_XS1_2*SIGe_1^2;
nXI3min(194,1) = 2*E_XF1_XS1_3*SIGe_1^2;
nXI3min(195,1) = 2*E_XF3_1*SIGe_1^2;
nXI3min(196,1) = 2*E_XF3_2*SIGe_1^2;
nXI3min(197,1) = 2*E_XF3_3*SIGe_1^2;
nXI3min(198,1) = 2*E_XF3_4*SIGe_1^2;
nXI3min(199,1) = 2*E_XF3_5*SIGe_1^2;
nXI3min(200,1) = 2*E_XF3_6*SIGe_1^2;
nXI3min(204,1) = 12*E_XF1_1*SIGe_1^3;
nXI3min(205,1) = 2*E_XF2_4*SIGe_1^2;
nXI3min(206,1) = 2*E_XF2_5*SIGe_1^2;
nXI3min(207,1) = 2*E_XF1_XS1_4*SIGe_1^2;
nXI3min(208,1) = 2*E_XF1_XS1_5*SIGe_1^2;
nXI3min(209,1) = 2*E_XF1_XS1_6*SIGe_1^2;
nXI3min(210,1) = 2*E_XF3_2*SIGe_1^2;
nXI3min(211,1) = 2*E_XF3_4*SIGe_1^2;
nXI3min(212,1) = 2*E_XF3_5*SIGe_1^2;
nXI3min(213,1) = 2*E_XF3_7*SIGe_1^2;
nXI3min(214,1) = 2*E_XF3_8*SIGe_1^2;
nXI3min(215,1) = 2*E_XF3_9*SIGe_1^2;
nXI3min(219,1) = 12*E_XF1_2*SIGe_1^3;
nXI3min(220,1) = 2*E_XF2_6*SIGe_1^2;
nXI3min(221,1) = 2*E_XF1_XS1_7*SIGe_1^2;
nXI3min(222,1) = 2*E_XF1_XS1_8*SIGe_1^2;
nXI3min(223,1) = 2*E_XF1_XS1_9*SIGe_1^2;
nXI3min(224,1) = 2*E_XF3_3*SIGe_1^2;
nXI3min(225,1) = 2*E_XF3_5*SIGe_1^2;
nXI3min(226,1) = 2*E_XF3_6*SIGe_1^2;
nXI3min(227,1) = 2*E_XF3_8*SIGe_1^2;
nXI3min(228,1) = 2*E_XF3_9*SIGe_1^2;
nXI3min(229,1) = 2*E_XF3_10*SIGe_1^2;
nXI3min(233,1) = 12*E_XF1_3*SIGe_1^3;
nXI3min(234,1) = 2*E_XS2_1*SIGe_1^2;
nXI3min(235,1) = 2*E_XS2_2*SIGe_1^2;
nXI3min(236,1) = 2*E_XS2_3*SIGe_1^2;
nXI3min(237,1) = 2*E_XF2_XS1_1*SIGe_1^2;
nXI3min(238,1) = 2*E_XF2_XS1_4*SIGe_1^2;
nXI3min(239,1) = 2*E_XF2_XS1_7*SIGe_1^2;
nXI3min(240,1) = 2*E_XF2_XS1_10*SIGe_1^2;
nXI3min(241,1) = 2*E_XF2_XS1_13*SIGe_1^2;
nXI3min(242,1) = 2*E_XF2_XS1_16*SIGe_1^2;
nXI3min(246,1) = 12*E_XS1_1*SIGe_1^3;
nXI3min(247,1) = 2*E_XS2_4*SIGe_1^2;
nXI3min(248,1) = 2*E_XS2_5*SIGe_1^2;
nXI3min(249,1) = 2*E_XF2_XS1_2*SIGe_1^2;
nXI3min(250,1) = 2*E_XF2_XS1_5*SIGe_1^2;
nXI3min(251,1) = 2*E_XF2_XS1_8*SIGe_1^2;
nXI3min(252,1) = 2*E_XF2_XS1_11*SIGe_1^2;
nXI3min(253,1) = 2*E_XF2_XS1_14*SIGe_1^2;
nXI3min(254,1) = 2*E_XF2_XS1_17*SIGe_1^2;
nXI3min(258,1) = 12*E_XS1_2*SIGe_1^3;
nXI3min(259,1) = 2*E_XS2_6*SIGe_1^2;
nXI3min(260,1) = 2*E_XF2_XS1_3*SIGe_1^2;
nXI3min(261,1) = 2*E_XF2_XS1_6*SIGe_1^2;
nXI3min(262,1) = 2*E_XF2_XS1_9*SIGe_1^2;
nXI3min(263,1) = 2*E_XF2_XS1_12*SIGe_1^2;
nXI3min(264,1) = 2*E_XF2_XS1_15*SIGe_1^2;
nXI3min(265,1) = 2*E_XF2_XS1_18*SIGe_1^2;
nXI3min(269,1) = 12*E_XS1_3*SIGe_1^3;
nXI3min(270,1) = 2*E_XF4_1*SIGe_1^2;
nXI3min(271,1) = 2*E_XF4_2*SIGe_1^2;
nXI3min(272,1) = 2*E_XF4_3*SIGe_1^2;
nXI3min(273,1) = 2*E_XF4_4*SIGe_1^2;
nXI3min(274,1) = 2*E_XF4_5*SIGe_1^2;
nXI3min(275,1) = 2*E_XF4_6*SIGe_1^2;
nXI3min(279,1) = 12*E_XF2_1*SIGe_1^3;
nXI3min(280,1) = 2*E_XF4_4*SIGe_1^2;
nXI3min(281,1) = 2*E_XF4_5*SIGe_1^2;
nXI3min(282,1) = 2*E_XF4_7*SIGe_1^2;
nXI3min(283,1) = 2*E_XF4_8*SIGe_1^2;
nXI3min(284,1) = 2*E_XF4_9*SIGe_1^2;
nXI3min(288,1) = 12*E_XF2_2*SIGe_1^3;
nXI3min(289,1) = 2*E_XF4_6*SIGe_1^2;
nXI3min(290,1) = 2*E_XF4_8*SIGe_1^2;
nXI3min(291,1) = 2*E_XF4_9*SIGe_1^2;
nXI3min(292,1) = 2*E_XF4_10*SIGe_1^2;
nXI3min(296,1) = 12*E_XF2_3*SIGe_1^3;
nXI3min(297,1) = 2*E_XF4_11*SIGe_1^2;
nXI3min(298,1) = 2*E_XF4_12*SIGe_1^2;
nXI3min(299,1) = 2*E_XF4_13*SIGe_1^2;
nXI3min(303,1) = 12*E_XF2_4*SIGe_1^3;
nXI3min(304,1) = 2*E_XF4_13*SIGe_1^2;
nXI3min(305,1) = 2*E_XF4_14*SIGe_1^2;
nXI3min(309,1) = 12*E_XF2_5*SIGe_1^3;
nXI3min(310,1) = 2*E_XF4_15*SIGe_1^2;
nXI3min(314,1) = 12*E_XF2_6*SIGe_1^3;
nXI3min(315,1) = 12*E_XF2_1*SIGe_1^3;
nXI3min(316,1) = 12*E_XF2_2*SIGe_1^3;
nXI3min(317,1) = 12*E_XF2_3*SIGe_1^3;
nXI3min(319,1) = 12*E_XF2_4*SIGe_1^3;
nXI3min(320,1) = 12*E_XF2_5*SIGe_1^3;
nXI3min(322,1) = 12*E_XF2_6*SIGe_1^3;
nXI3min(324,1) = 90*SIGe_1^4;
nXI3min(337,1) = 3*E_XF3_1*SIGe_1^2;
nXI3min(338,1) = 3*E_XF3_2*SIGe_1^2;
nXI3min(339,1) = 3*E_XF3_3*SIGe_1^2;
nXI3min(352,1) = 3*E_XF3_2*SIGe_1^2;
nXI3min(353,1) = 3*E_XF3_4*SIGe_1^2;
nXI3min(354,1) = 3*E_XF3_5*SIGe_1^2;
nXI3min(366,1) = 3*E_XF3_3*SIGe_1^2;
nXI3min(367,1) = 3*E_XF3_5*SIGe_1^2;
nXI3min(368,1) = 3*E_XF3_6*SIGe_1^2;
nXI3min(379,1) = 3*E_XF2_XS1_1*SIGe_1^2;
nXI3min(380,1) = 3*E_XF2_XS1_4*SIGe_1^2;
nXI3min(381,1) = 3*E_XF2_XS1_7*SIGe_1^2;
nXI3min(391,1) = 3*E_XF2_XS1_2*SIGe_1^2;
nXI3min(392,1) = 3*E_XF2_XS1_5*SIGe_1^2;
nXI3min(393,1) = 3*E_XF2_XS1_8*SIGe_1^2;
nXI3min(402,1) = 3*E_XF2_XS1_3*SIGe_1^2;
nXI3min(403,1) = 3*E_XF2_XS1_6*SIGe_1^2;
nXI3min(404,1) = 3*E_XF2_XS1_9*SIGe_1^2;
nXI3min(412,1) = 3*E_XF4_1*SIGe_1^2;
nXI3min(413,1) = 3*E_XF4_2*SIGe_1^2;
nXI3min(414,1) = 3*E_XF4_3*SIGe_1^2;
nXI3min(421,1) = 3*E_XF4_2*SIGe_1^2;
nXI3min(422,1) = 3*E_XF4_4*SIGe_1^2;
nXI3min(423,1) = 3*E_XF4_5*SIGe_1^2;
nXI3min(429,1) = 3*E_XF4_3*SIGe_1^2;
nXI3min(430,1) = 3*E_XF4_5*SIGe_1^2;
nXI3min(431,1) = 3*E_XF4_6*SIGe_1^2;
nXI3min(436,1) = 3*E_XF4_4*SIGe_1^2;
nXI3min(437,1) = 3*E_XF4_7*SIGe_1^2;
nXI3min(438,1) = 3*E_XF4_8*SIGe_1^2;
nXI3min(442,1) = 3*E_XF4_5*SIGe_1^2;
nXI3min(443,1) = 3*E_XF4_8*SIGe_1^2;
nXI3min(444,1) = 3*E_XF4_9*SIGe_1^2;
nXI3min(447,1) = 3*E_XF4_6*SIGe_1^2;
nXI3min(448,1) = 3*E_XF4_9*SIGe_1^2;
nXI3min(449,1) = 3*E_XF4_10*SIGe_1^2;
nXI3min(454,1) = 15*E_XF2_1*SIGe_1^3;
nXI3min(457,1) = 15*E_XF2_2*SIGe_1^3;
nXI3min(459,1) = 15*E_XF2_3*SIGe_1^3;
nXI3min(472,1) = 3*E_XF3_4*SIGe_1^2;
nXI3min(473,1) = 3*E_XF3_7*SIGe_1^2;
nXI3min(474,1) = 3*E_XF3_8*SIGe_1^2;
nXI3min(486,1) = 3*E_XF3_5*SIGe_1^2;
nXI3min(487,1) = 3*E_XF3_8*SIGe_1^2;
nXI3min(488,1) = 3*E_XF3_9*SIGe_1^2;
nXI3min(499,1) = 3*E_XF2_XS1_4*SIGe_1^2;
nXI3min(500,1) = 3*E_XF2_XS1_10*SIGe_1^2;
nXI3min(501,1) = 3*E_XF2_XS1_13*SIGe_1^2;
nXI3min(511,1) = 3*E_XF2_XS1_5*SIGe_1^2;
nXI3min(512,1) = 3*E_XF2_XS1_11*SIGe_1^2;
nXI3min(513,1) = 3*E_XF2_XS1_14*SIGe_1^2;
nXI3min(522,1) = 3*E_XF2_XS1_6*SIGe_1^2;
nXI3min(523,1) = 3*E_XF2_XS1_12*SIGe_1^2;
nXI3min(524,1) = 3*E_XF2_XS1_15*SIGe_1^2;
nXI3min(532,1) = 3*E_XF4_2*SIGe_1^2;
nXI3min(533,1) = 3*E_XF4_4*SIGe_1^2;
nXI3min(534,1) = 3*E_XF4_5*SIGe_1^2;
nXI3min(541,1) = 3*E_XF4_4*SIGe_1^2;
nXI3min(542,1) = 3*E_XF4_7*SIGe_1^2;
nXI3min(543,1) = 3*E_XF4_8*SIGe_1^2;
nXI3min(549,1) = 3*E_XF4_5*SIGe_1^2;
nXI3min(550,1) = 3*E_XF4_8*SIGe_1^2;
nXI3min(551,1) = 3*E_XF4_9*SIGe_1^2;
nXI3min(556,1) = 3*E_XF4_7*SIGe_1^2;
nXI3min(557,1) = 3*E_XF4_11*SIGe_1^2;
nXI3min(558,1) = 3*E_XF4_12*SIGe_1^2;
nXI3min(562,1) = 3*E_XF4_8*SIGe_1^2;
nXI3min(563,1) = 3*E_XF4_12*SIGe_1^2;
nXI3min(564,1) = 3*E_XF4_13*SIGe_1^2;
nXI3min(567,1) = 3*E_XF4_9*SIGe_1^2;
nXI3min(568,1) = 3*E_XF4_13*SIGe_1^2;
nXI3min(569,1) = 3*E_XF4_14*SIGe_1^2;
nXI3min(574,1) = 15*E_XF2_2*SIGe_1^3;
nXI3min(577,1) = 15*E_XF2_4*SIGe_1^3;
nXI3min(579,1) = 15*E_XF2_5*SIGe_1^3;
nXI3min(591,1) = 3*E_XF3_6*SIGe_1^2;
nXI3min(592,1) = 3*E_XF3_9*SIGe_1^2;
nXI3min(593,1) = 3*E_XF3_10*SIGe_1^2;
nXI3min(604,1) = 3*E_XF2_XS1_7*SIGe_1^2;
nXI3min(605,1) = 3*E_XF2_XS1_13*SIGe_1^2;
nXI3min(606,1) = 3*E_XF2_XS1_16*SIGe_1^2;
nXI3min(616,1) = 3*E_XF2_XS1_8*SIGe_1^2;
nXI3min(617,1) = 3*E_XF2_XS1_14*SIGe_1^2;
nXI3min(618,1) = 3*E_XF2_XS1_17*SIGe_1^2;
nXI3min(627,1) = 3*E_XF2_XS1_9*SIGe_1^2;
nXI3min(628,1) = 3*E_XF2_XS1_15*SIGe_1^2;
nXI3min(629,1) = 3*E_XF2_XS1_18*SIGe_1^2;
nXI3min(637,1) = 3*E_XF4_3*SIGe_1^2;
nXI3min(638,1) = 3*E_XF4_5*SIGe_1^2;
nXI3min(639,1) = 3*E_XF4_6*SIGe_1^2;
nXI3min(646,1) = 3*E_XF4_5*SIGe_1^2;
nXI3min(647,1) = 3*E_XF4_8*SIGe_1^2;
nXI3min(648,1) = 3*E_XF4_9*SIGe_1^2;
nXI3min(654,1) = 3*E_XF4_6*SIGe_1^2;
nXI3min(655,1) = 3*E_XF4_9*SIGe_1^2;
nXI3min(656,1) = 3*E_XF4_10*SIGe_1^2;
nXI3min(661,1) = 3*E_XF4_8*SIGe_1^2;
nXI3min(662,1) = 3*E_XF4_12*SIGe_1^2;
nXI3min(663,1) = 3*E_XF4_13*SIGe_1^2;
nXI3min(667,1) = 3*E_XF4_9*SIGe_1^2;
nXI3min(668,1) = 3*E_XF4_13*SIGe_1^2;
nXI3min(669,1) = 3*E_XF4_14*SIGe_1^2;
nXI3min(672,1) = 3*E_XF4_10*SIGe_1^2;
nXI3min(673,1) = 3*E_XF4_14*SIGe_1^2;
nXI3min(674,1) = 3*E_XF4_15*SIGe_1^2;
nXI3min(679,1) = 15*E_XF2_3*SIGe_1^3;
nXI3min(682,1) = 15*E_XF2_5*SIGe_1^3;
nXI3min(684,1) = 15*E_XF2_6*SIGe_1^3;
nXI3min(695,1) = 3*E_XF1_XS2_1*SIGe_1^2;
nXI3min(696,1) = 3*E_XF1_XS2_7*SIGe_1^2;
nXI3min(697,1) = 3*E_XF1_XS2_13*SIGe_1^2;
nXI3min(707,1) = 3*E_XF1_XS2_2*SIGe_1^2;
nXI3min(708,1) = 3*E_XF1_XS2_8*SIGe_1^2;
nXI3min(709,1) = 3*E_XF1_XS2_14*SIGe_1^2;
nXI3min(718,1) = 3*E_XF1_XS2_3*SIGe_1^2;
nXI3min(719,1) = 3*E_XF1_XS2_9*SIGe_1^2;
nXI3min(720,1) = 3*E_XF1_XS2_15*SIGe_1^2;
nXI3min(728,1) = 3*E_XF3_XS1_1*SIGe_1^2;
nXI3min(729,1) = 3*E_XF3_XS1_4*SIGe_1^2;
nXI3min(730,1) = 3*E_XF3_XS1_7*SIGe_1^2;
nXI3min(737,1) = 3*E_XF3_XS1_4*SIGe_1^2;
nXI3min(738,1) = 3*E_XF3_XS1_10*SIGe_1^2;
nXI3min(739,1) = 3*E_XF3_XS1_13*SIGe_1^2;
nXI3min(745,1) = 3*E_XF3_XS1_7*SIGe_1^2;
nXI3min(746,1) = 3*E_XF3_XS1_13*SIGe_1^2;
nXI3min(747,1) = 3*E_XF3_XS1_16*SIGe_1^2;
nXI3min(752,1) = 3*E_XF3_XS1_10*SIGe_1^2;
nXI3min(753,1) = 3*E_XF3_XS1_19*SIGe_1^2;
nXI3min(754,1) = 3*E_XF3_XS1_22*SIGe_1^2;
nXI3min(758,1) = 3*E_XF3_XS1_13*SIGe_1^2;
nXI3min(759,1) = 3*E_XF3_XS1_22*SIGe_1^2;
nXI3min(760,1) = 3*E_XF3_XS1_25*SIGe_1^2;
nXI3min(763,1) = 3*E_XF3_XS1_16*SIGe_1^2;
nXI3min(764,1) = 3*E_XF3_XS1_25*SIGe_1^2;
nXI3min(765,1) = 3*E_XF3_XS1_28*SIGe_1^2;
nXI3min(770,1) = 15*E_XF1_XS1_1*SIGe_1^3;
nXI3min(773,1) = 15*E_XF1_XS1_4*SIGe_1^3;
nXI3min(775,1) = 15*E_XF1_XS1_7*SIGe_1^3;
nXI3min(785,1) = 3*E_XF1_XS2_4*SIGe_1^2;
nXI3min(786,1) = 3*E_XF1_XS2_10*SIGe_1^2;
nXI3min(787,1) = 3*E_XF1_XS2_16*SIGe_1^2;
nXI3min(796,1) = 3*E_XF1_XS2_5*SIGe_1^2;
nXI3min(797,1) = 3*E_XF1_XS2_11*SIGe_1^2;
nXI3min(798,1) = 3*E_XF1_XS2_17*SIGe_1^2;
nXI3min(806,1) = 3*E_XF3_XS1_2*SIGe_1^2;
nXI3min(807,1) = 3*E_XF3_XS1_5*SIGe_1^2;
nXI3min(808,1) = 3*E_XF3_XS1_8*SIGe_1^2;
nXI3min(815,1) = 3*E_XF3_XS1_5*SIGe_1^2;
nXI3min(816,1) = 3*E_XF3_XS1_11*SIGe_1^2;
nXI3min(817,1) = 3*E_XF3_XS1_14*SIGe_1^2;
nXI3min(823,1) = 3*E_XF3_XS1_8*SIGe_1^2;
nXI3min(824,1) = 3*E_XF3_XS1_14*SIGe_1^2;
nXI3min(825,1) = 3*E_XF3_XS1_17*SIGe_1^2;
nXI3min(830,1) = 3*E_XF3_XS1_11*SIGe_1^2;
nXI3min(831,1) = 3*E_XF3_XS1_20*SIGe_1^2;
nXI3min(832,1) = 3*E_XF3_XS1_23*SIGe_1^2;
nXI3min(836,1) = 3*E_XF3_XS1_14*SIGe_1^2;
nXI3min(837,1) = 3*E_XF3_XS1_23*SIGe_1^2;
nXI3min(838,1) = 3*E_XF3_XS1_26*SIGe_1^2;
nXI3min(841,1) = 3*E_XF3_XS1_17*SIGe_1^2;
nXI3min(842,1) = 3*E_XF3_XS1_26*SIGe_1^2;
nXI3min(843,1) = 3*E_XF3_XS1_29*SIGe_1^2;
nXI3min(848,1) = 15*E_XF1_XS1_2*SIGe_1^3;
nXI3min(851,1) = 15*E_XF1_XS1_5*SIGe_1^3;
nXI3min(853,1) = 15*E_XF1_XS1_8*SIGe_1^3;
nXI3min(862,1) = 3*E_XF1_XS2_6*SIGe_1^2;
nXI3min(863,1) = 3*E_XF1_XS2_12*SIGe_1^2;
nXI3min(864,1) = 3*E_XF1_XS2_18*SIGe_1^2;
nXI3min(872,1) = 3*E_XF3_XS1_3*SIGe_1^2;
nXI3min(873,1) = 3*E_XF3_XS1_6*SIGe_1^2;
nXI3min(874,1) = 3*E_XF3_XS1_9*SIGe_1^2;
nXI3min(881,1) = 3*E_XF3_XS1_6*SIGe_1^2;
nXI3min(882,1) = 3*E_XF3_XS1_12*SIGe_1^2;
nXI3min(883,1) = 3*E_XF3_XS1_15*SIGe_1^2;
nXI3min(889,1) = 3*E_XF3_XS1_9*SIGe_1^2;
nXI3min(890,1) = 3*E_XF3_XS1_15*SIGe_1^2;
nXI3min(891,1) = 3*E_XF3_XS1_18*SIGe_1^2;
nXI3min(896,1) = 3*E_XF3_XS1_12*SIGe_1^2;
nXI3min(897,1) = 3*E_XF3_XS1_21*SIGe_1^2;
nXI3min(898,1) = 3*E_XF3_XS1_24*SIGe_1^2;
nXI3min(902,1) = 3*E_XF3_XS1_15*SIGe_1^2;
nXI3min(903,1) = 3*E_XF3_XS1_24*SIGe_1^2;
nXI3min(904,1) = 3*E_XF3_XS1_27*SIGe_1^2;
nXI3min(907,1) = 3*E_XF3_XS1_18*SIGe_1^2;
nXI3min(908,1) = 3*E_XF3_XS1_27*SIGe_1^2;
nXI3min(909,1) = 3*E_XF3_XS1_30*SIGe_1^2;
nXI3min(914,1) = 15*E_XF1_XS1_3*SIGe_1^3;
nXI3min(917,1) = 15*E_XF1_XS1_6*SIGe_1^3;
nXI3min(919,1) = 15*E_XF1_XS1_9*SIGe_1^3;
nXI3min(927,1) = 3*E_XF5_1*SIGe_1^2;
nXI3min(928,1) = 3*E_XF5_2*SIGe_1^2;
nXI3min(929,1) = 3*E_XF5_3*SIGe_1^2;
nXI3min(936,1) = 3*E_XF5_2*SIGe_1^2;
nXI3min(937,1) = 3*E_XF5_4*SIGe_1^2;
nXI3min(938,1) = 3*E_XF5_5*SIGe_1^2;
nXI3min(944,1) = 3*E_XF5_3*SIGe_1^2;
nXI3min(945,1) = 3*E_XF5_5*SIGe_1^2;
nXI3min(946,1) = 3*E_XF5_6*SIGe_1^2;
nXI3min(951,1) = 3*E_XF5_4*SIGe_1^2;
nXI3min(952,1) = 3*E_XF5_7*SIGe_1^2;
nXI3min(953,1) = 3*E_XF5_8*SIGe_1^2;
nXI3min(957,1) = 3*E_XF5_5*SIGe_1^2;
nXI3min(958,1) = 3*E_XF5_8*SIGe_1^2;
nXI3min(959,1) = 3*E_XF5_9*SIGe_1^2;
nXI3min(962,1) = 3*E_XF5_6*SIGe_1^2;
nXI3min(963,1) = 3*E_XF5_9*SIGe_1^2;
nXI3min(964,1) = 3*E_XF5_10*SIGe_1^2;
nXI3min(969,1) = 15*E_XF3_1*SIGe_1^3;
nXI3min(972,1) = 15*E_XF3_2*SIGe_1^3;
nXI3min(974,1) = 15*E_XF3_3*SIGe_1^3;
nXI3min(981,1) = 3*E_XF5_4*SIGe_1^2;
nXI3min(982,1) = 3*E_XF5_7*SIGe_1^2;
nXI3min(983,1) = 3*E_XF5_8*SIGe_1^2;
nXI3min(989,1) = 3*E_XF5_5*SIGe_1^2;
nXI3min(990,1) = 3*E_XF5_8*SIGe_1^2;
nXI3min(991,1) = 3*E_XF5_9*SIGe_1^2;
nXI3min(996,1) = 3*E_XF5_7*SIGe_1^2;
nXI3min(997,1) = 3*E_XF5_11*SIGe_1^2;
nXI3min(998,1) = 3*E_XF5_12*SIGe_1^2;
nXI3min(1002,1) = 3*E_XF5_8*SIGe_1^2;
nXI3min(1003,1) = 3*E_XF5_12*SIGe_1^2;
nXI3min(1004,1) = 3*E_XF5_13*SIGe_1^2;
nXI3min(1007,1) = 3*E_XF5_9*SIGe_1^2;
nXI3min(1008,1) = 3*E_XF5_13*SIGe_1^2;
nXI3min(1009,1) = 3*E_XF5_14*SIGe_1^2;
nXI3min(1014,1) = 15*E_XF3_2*SIGe_1^3;
nXI3min(1017,1) = 15*E_XF3_4*SIGe_1^3;
nXI3min(1019,1) = 15*E_XF3_5*SIGe_1^3;
nXI3min(1025,1) = 3*E_XF5_6*SIGe_1^2;
nXI3min(1026,1) = 3*E_XF5_9*SIGe_1^2;
nXI3min(1027,1) = 3*E_XF5_10*SIGe_1^2;
nXI3min(1032,1) = 3*E_XF5_8*SIGe_1^2;
nXI3min(1033,1) = 3*E_XF5_12*SIGe_1^2;
nXI3min(1034,1) = 3*E_XF5_13*SIGe_1^2;
nXI3min(1038,1) = 3*E_XF5_9*SIGe_1^2;
nXI3min(1039,1) = 3*E_XF5_13*SIGe_1^2;
nXI3min(1040,1) = 3*E_XF5_14*SIGe_1^2;
nXI3min(1043,1) = 3*E_XF5_10*SIGe_1^2;
nXI3min(1044,1) = 3*E_XF5_14*SIGe_1^2;
nXI3min(1045,1) = 3*E_XF5_15*SIGe_1^2;
nXI3min(1050,1) = 15*E_XF3_3*SIGe_1^3;
nXI3min(1053,1) = 15*E_XF3_5*SIGe_1^3;
nXI3min(1055,1) = 15*E_XF3_6*SIGe_1^3;
nXI3min(1060,1) = 3*E_XF5_11*SIGe_1^2;
nXI3min(1061,1) = 3*E_XF5_16*SIGe_1^2;
nXI3min(1062,1) = 3*E_XF5_17*SIGe_1^2;
nXI3min(1066,1) = 3*E_XF5_12*SIGe_1^2;
nXI3min(1067,1) = 3*E_XF5_17*SIGe_1^2;
nXI3min(1068,1) = 3*E_XF5_18*SIGe_1^2;
nXI3min(1071,1) = 3*E_XF5_13*SIGe_1^2;
nXI3min(1072,1) = 3*E_XF5_18*SIGe_1^2;
nXI3min(1073,1) = 3*E_XF5_19*SIGe_1^2;
nXI3min(1078,1) = 15*E_XF3_4*SIGe_1^3;
nXI3min(1081,1) = 15*E_XF3_7*SIGe_1^3;
nXI3min(1083,1) = 15*E_XF3_8*SIGe_1^3;
nXI3min(1087,1) = 3*E_XF5_13*SIGe_1^2;
nXI3min(1088,1) = 3*E_XF5_18*SIGe_1^2;
nXI3min(1089,1) = 3*E_XF5_19*SIGe_1^2;
nXI3min(1092,1) = 3*E_XF5_14*SIGe_1^2;
nXI3min(1093,1) = 3*E_XF5_19*SIGe_1^2;
nXI3min(1094,1) = 3*E_XF5_20*SIGe_1^2;
nXI3min(1099,1) = 15*E_XF3_5*SIGe_1^3;
nXI3min(1102,1) = 15*E_XF3_8*SIGe_1^3;
nXI3min(1104,1) = 15*E_XF3_9*SIGe_1^3;
nXI3min(1107,1) = 3*E_XF5_15*SIGe_1^2;
nXI3min(1108,1) = 3*E_XF5_20*SIGe_1^2;
nXI3min(1109,1) = 3*E_XF5_21*SIGe_1^2;
nXI3min(1114,1) = 15*E_XF3_6*SIGe_1^3;
nXI3min(1117,1) = 15*E_XF3_9*SIGe_1^3;
nXI3min(1119,1) = 15*E_XF3_10*SIGe_1^3;
nXI3min(1121,1) = 15*E_XF3_1*SIGe_1^3;
nXI3min(1122,1) = 15*E_XF3_2*SIGe_1^3;
nXI3min(1123,1) = 15*E_XF3_3*SIGe_1^3;
nXI3min(1125,1) = 15*E_XF3_4*SIGe_1^3;
nXI3min(1126,1) = 15*E_XF3_5*SIGe_1^3;
nXI3min(1128,1) = 15*E_XF3_6*SIGe_1^3;
nXI3min(1130,1) = 105*E_XF1_1*SIGe_1^4;
nXI3min(1131,1) = 15*E_XF3_7*SIGe_1^3;
nXI3min(1132,1) = 15*E_XF3_8*SIGe_1^3;
nXI3min(1134,1) = 15*E_XF3_9*SIGe_1^3;
nXI3min(1136,1) = 105*E_XF1_2*SIGe_1^4;
nXI3min(1137,1) = 15*E_XF3_10*SIGe_1^3;
nXI3min(1139,1) = 105*E_XF1_3*SIGe_1^4;
