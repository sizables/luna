lune = window.lune = require 'lune'
window.onload = ->
  console.log lune.phase()

# Phaser = require 'phaser'
# window.onload = ->
#   game = new Phaser.Game 800, 600, Phaser.AUTO, '', {
#     preload: ->
#       game.load.image 'logo', 'luna.png'
#     create: ->
#       logo = game.add.sprite game.world.centerX, game.world.centerY, 'logo'
#       logo.anchor.setTo 0.5, 0.5
#   }
