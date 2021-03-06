#!/usr/bin/env python
#***************************************************************************
#
#   Copyright (c) 2015 PX4 Development Team. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
# 3. Neither the name PX4 nor the names of its contributors may be
#    used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
#***************************************************************************/

#
# @author Andreas Antener <andreas@uaventure.com>
#
PKG = 'px4'

import unittest
import rospy
import math
import rosbag
import sys
import os
import time

import mavros
from pymavlink import mavutil
from mavros import mavlink

from geometry_msgs.msg import PoseStamped
from mavros_msgs.srv import CommandLong, WaypointPush
from mavros_msgs.msg import Mavlink, Waypoint, ExtendedState
from sensor_msgs.msg import NavSatFix
from mavros.mission import QGroundControlWP
#from px4_test_helper import PX4TestHelper

class MavrosMissionTest(unittest.TestCase):
    """
    Run a mission
    """

    def setUp(self):
        rospy.init_node('test_node', anonymous=True)

        self.rate = rospy.Rate(10) # 10hz
        self.has_global_pos = False
        self.local_position = PoseStamped()
        self.global_position = NavSatFix()
        self.extended_state = ExtendedState()
        self.home_alt = 0
        self.mc_rad = 5
        self.fw_rad = 60
        self.fw_alt_rad = 10
        self.last_alt_d = 9999
        self.last_pos_d = 9999
        self.mission_name = ""

        # need to simulate heartbeat for datalink loss detection
        rospy.Timer(rospy.Duration(0.5), self.send_heartbeat)

        rospy.wait_for_service('mavros/cmd/command', 30)
        self.pub_mavlink = rospy.Publisher('mavlink/to', Mavlink, queue_size=1)
        self._srv_cmd_long = rospy.ServiceProxy('mavros/cmd/command', CommandLong, persistent=True)
        self._srv_wp_push = rospy.ServiceProxy('mavros/mission/push', WaypointPush, persistent=True)

        rospy.Subscriber("mavros/local_position/pose", PoseStamped, self.position_callback)
        rospy.Subscriber("mavros/global_position/global", NavSatFix, self.global_position_callback)
        rospy.Subscriber("mavros/extended_state", ExtendedState, self.extended_state_callback)

    def tearDown(self):
        #self.helper.tearDown()
        pass

    #
    # General callback functions used in tests
    #
    def position_callback(self, data):
        self.local_position = data

    def global_position_callback(self, data):
        self.global_position = data

        if not self.has_global_pos:
            self.home_alt = data.altitude
            self.has_global_pos = True

    def extended_state_callback(self, data):
        self.extended_state = data

    #
    # Helper methods
    #
    def is_at_position(self, lat, lon, alt, xy_offset, z_offset):
        R = 6371000 # metres
        rlat1 = math.radians(lat)
        rlat2 = math.radians(self.global_position.latitude)

        rlat_d = math.radians(self.global_position.latitude - lat)
        rlon_d = math.radians(self.global_position.longitude - lon)

        a = (math.sin(rlat_d / 2) * math.sin(rlat_d / 2) +
             math.cos(rlat1) * math.cos(rlat2) *
             math.sin(rlon_d / 2) * math.sin(rlon_d / 2))
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

        d = R * c
        alt_d = abs(alt - self.global_position.altitude)

        #rospy.loginfo("d: %f, alt_d: %f", d, alt_d)

        # remember best distances
        if self.last_pos_d > d:
            self.last_pos_d = d
        if self.last_alt_d > alt_d:
            self.last_alt_d = alt_d

        return d < xy_offset and alt_d < z_offset

    def reach_position(self, lat, lon, alt, timeout, index):
        # reset best distances
        self.last_alt_d = 9999
        self.last_pos_d = 9999

        rospy.loginfo("trying to reach waypoint " +
            "lat: %f, lon: %f, alt: %f, timeout: %d, index: %d" %
            (lat, lon, alt, timeout, index))

        # does it reach the position in X seconds?
        count = 0
        while count < timeout:
            # use MC radius by default
            # FIXME: also check MAV_TYPE from system status, otherwise pure fixed-wing won't work
            xy_radius = self.mc_rad
            z_radius = self.mc_rad

            # use FW radius if in FW or in transition
            if (self.extended_state.vtol_state == ExtendedState.VTOL_STATE_FW or
                    self.extended_state.vtol_state == ExtendedState.VTOL_STATE_TRANSITION_TO_MC or
                    self.extended_state.vtol_state == ExtendedState.VTOL_STATE_TRANSITION_TO_FW):
                xy_radius = self.fw_rad
                z_radius = self.fw_alt_rad

            if self.is_at_position(lat, lon, alt, xy_radius, z_radius):
                rospy.loginfo("position reached, index: %d, count: %d, pos_d: %f, alt_d: %f" %
                    (index, count, self.last_pos_d, self.last_alt_d))
                break

            count = count + 1
            self.rate.sleep()

        self.assertTrue(count < timeout, (("(%s) took too long to get to position " +
            "lat: %f, lon: %f, alt: %f, xy off: %f, z off: %f, timeout: %d, index: %d, pos_d: %f, alt_d: %f") %
            (self.mission_name, lat, lon, alt, xy_radius, z_radius, timeout, index, self.last_pos_d, self.last_alt_d)))

    def run_mission(self):
        """switch mode: auto and arm"""
        self._srv_cmd_long(False, 176, False,
                           # custom, auto, mission
                           1, 4, 4, 0, 0, 0, 0)
        # make sure the first command doesn't get lost
        time.sleep(1)

        self._srv_cmd_long(False, 400, False,
                           # arm
                           1, 0, 0, 0, 0, 0, 0)

    def wait_until_ready(self):
        """FIXME: hack to wait for simulation to be ready"""
        while not self.has_global_pos:
            self.rate.sleep()

    def wait_on_landing(self, timeout, index):
        """Wait for landed state"""

        rospy.loginfo("waiting for landing " +
            "timeout: %d, index: %d" %
            (timeout, index))

        count = 0
        while count < timeout:
            if self.extended_state.landed_state == ExtendedState.LANDED_STATE_ON_GROUND:
                break

            count = count + 1
            self.rate.sleep()

        self.assertTrue(count < timeout, ("(%s) landing not detected after landing WP " +
            "timeout: %d, index: %d") %
            (self.mission_name, timeout, index))

    def wait_on_transition(self, transition, timeout, index):
        """Wait for VTOL transition"""

        rospy.loginfo("waiting for VTOL transition " +
            "transition: %d, timeout: %d, index: %d" %
            (transition, timeout, index))

        count = 0
        while count < timeout:
            # transition to MC
            if (transition == ExtendedState.VTOL_STATE_MC and
                    self.extended_state.vtol_state == ExtendedState.VTOL_STATE_MC):
                break

            # transition to FW
            if (transition == ExtendedState.VTOL_STATE_FW and
                    self.extended_state.vtol_state == ExtendedState.VTOL_STATE_FW):
                break

            count = count + 1
            self.rate.sleep()

        self.assertTrue(count < timeout, ("(%s) transition not detected " +
            "timeout: %d, index: %d") %
            (self.mission_name, timeout, index))

    def send_heartbeat(self, event=None):
        # mav type gcs
        mavmsg = mavutil.mavlink.MAVLink_heartbeat_message(6, 0, 0, 0, 0, 0)
        # XXX: hack: using header object to set mav properties
        mavmsg.pack(mavutil.mavlink.MAVLink_header(0, 0, 0, 2, 1))
        rosmsg = mavlink.convert_to_rosmsg(mavmsg)
        self.pub_mavlink.publish(rosmsg)

    def test_mission(self):
        """Test mission"""

        if len(sys.argv) < 2:
            self.fail("usage: mission_test.py mission_file")
            return

        self.mission_name = sys.argv[1]
        mission_file = os.path.dirname(os.path.realpath(__file__)) + "/" + sys.argv[1]

        rospy.loginfo("reading mission %s", mission_file)
        mission = QGroundControlWP()
        wps = []
        for waypoint in mission.read(open(mission_file, 'r')):
            wps.append(waypoint)
            rospy.logdebug(waypoint)

        rospy.loginfo("wait until ready")
        self.wait_until_ready()

        rospy.loginfo("send mission")
        res = self._srv_wp_push(wps)
        rospy.loginfo(res)
        self.assertTrue(res.success, "(%s) mission could not be transfered" % self.mission_name)

        rospy.loginfo("run mission")
        self.run_mission()

        index = 0
        for waypoint in wps:
            # only check position for waypoints where this makes sense
            if waypoint.frame == Waypoint.FRAME_GLOBAL_REL_ALT or waypoint.frame == Waypoint.FRAME_GLOBAL:
                alt = waypoint.z_alt
                if waypoint.frame == Waypoint.FRAME_GLOBAL_REL_ALT:
                    alt += self.home_alt
                self.reach_position(waypoint.x_lat, waypoint.y_long, alt, 600, index)

            # check if VTOL transition happens if applicable
            if waypoint.command == 84 or waypoint.command == 85 or waypoint.command == 3000:
                transition = waypoint.param1

                if waypoint.command == 84: # VTOL takeoff implies transition to FW
                    transition = ExtendedState.VTOL_STATE_FW

                if waypoint.command == 85: # VTOL takeoff implies transition to MC
                    transition = ExtendedState.VTOL_STATE_MC

                self.wait_on_transition(transition, 600, index)

            # after reaching position, wait for landing detection if applicable
            if waypoint.command == 85 or waypoint.command == 21:
                self.wait_on_landing(600, index)

            index += 1


if __name__ == '__main__':
    import rostest
    name = "mavros_mission_test"
    if len(sys.argv) > 1:
        name += "-%s" % sys.argv[1]
    rostest.rosrun(PKG, name, MavrosMissionTest)
