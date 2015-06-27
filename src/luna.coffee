lune = require 'lune'
console.log lune.phase()

require './gif.coffee'

class Luna extends Phaser.Game
  constructor: ->
    unless @ instanceof Luna then return new Luna
    @game = super window.innerWidth, window.innerHeight, Phaser.AUTO, '', @
  preload: ->
    @game.load.crossOrigin = '*'
    @game.load.spritegif 'm8hYEjb', 'http://i.imgur.com/m8hYEjb.gif'

  create: ->

window.game = new Luna
